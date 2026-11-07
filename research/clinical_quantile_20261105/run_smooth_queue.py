import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import pandas as pd
import yaml

config_path = Path(sys.argv[1] if len(sys.argv) > 1 else "config/clinical_quantile_smooth_20261106.yml")
config = yaml.safe_load(config_path.read_text())
root = Path(config["output_root"])
if not root.exists():
    root.mkdir(parents=True)
    implementation = root / "implementation"
    implementation.mkdir()
    for source in [config_path, Path("research/clinical_quantile_20261105/engine.R"),
                   Path("research/clinical_quantile_20261105/run_smooth_fold.R"),
                   Path("research/clinical_quantile_20261105/smooth_path.py"),
                   Path("research/clinical_search_20261105/panel.R"),
                   Path("research/empirical_selection_20261101/prediction_core.R")]:
        shutil.copy2(source, implementation / source.name)
    for path in config["source_studies"]:
        source = Path(path)
        shutil.copy2(source, implementation / source.name)
    (root / "request.json").write_text(json.dumps({
        "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "method_id": config["method_id"], "formal_validation_used": False}, indent=2) + "\n")
implementation = root / "implementation"
assert yaml.safe_load((implementation / config_path.name).read_text()) == config
source_paths = [implementation / Path(path).name for path in config["source_studies"]]
studies = [yaml.safe_load(path.read_text()) for path in source_paths]
pending = []
for fold in range(1, max(study["outer_folds"] for study in studies) + 1):
    for path, study in zip(source_paths, studies):
        if fold > study["outer_folds"]:
            continue
        directory = root / study["study_id"]
        directory.mkdir(exist_ok=True)
        missing = [tau for tau in config["quantiles"] if not
                   (directory / f"outer_{fold}" / f"tau{round(tau * 100):03d}" / "status.json").exists()]
        if missing:
            pending.append((path, study["study_id"], fold))
environment = os.environ.copy()
environment.update({"TMPDIR": str(Path("results/scratch/runtime").resolve()),
    "RETICULATE_PYTHON": str(Path.cwd() / "results/scratch/cq_solver_env/bin/python"),
    "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1",
    "NUMBA_NUM_THREADS": "1", "R_MAX_VSIZE": "24Gb"})
running = []
while True:
    while pending and len(running) < config["workers"]:
        path, study, fold = pending.pop(0)
        log = (root / f"{study}_fold_{fold}.log").open("a")
        process = subprocess.Popen(["Rscript", str(implementation / "run_smooth_fold.R"),
                                    str(implementation / config_path.name), str(path), str(fold)],
                                    env=environment, stdout=log, stderr=subprocess.STDOUT)
        running.append((process, log, study, fold))
        print(f"开始计算 {study} fold {fold} PID {process.pid}", flush=True)
    remaining = []
    for process, log, study, fold in running:
        code = process.poll()
        if code is None:
            remaining.append((process, log, study, fold))
            continue
        log.close()
        if code != 0:
            for other, stream, _, _ in running:
                if other.poll() is None:
                    other.terminate()
                stream.close()
            (root / "status.json").write_text(json.dumps({"status": "failed", "study_id": study,
                "fold": fold, "exit_code": code}, indent=2) + "\n")
            raise RuntimeError(f"{study} fold {fold} 计算终止，退出码 {code}")
        print(f"完成计算 {study} fold {fold}", flush=True)
    running = remaining
    progress = []
    for study in studies:
        directory = root / study["study_id"]
        files = sorted(directory.glob("outer_*/tau*/status.json"))
        for file in files:
            assert json.loads(file.read_text())["status"] == "completed"
        complete = len(files) == study["outer_folds"] * len(config["quantiles"])
        state = {"status": "completed" if complete else "running", "fits": len(files),
                 "formal_validation_used": False, "manuscript_adoption": False}
        (directory / "status.json").write_text(json.dumps(state, indent=2) + "\n")
        if complete and not (directory / "predictions.csv").exists():
            pd.concat([pd.read_csv(file.parent / "prediction.csv") for file in files], ignore_index=True).to_csv(directory / "predictions.csv", index=False)
            pd.concat([pd.read_csv(file.parent / "diagnostics.csv") for file in files], ignore_index=True).to_csv(directory / "diagnostics.csv", index=False)
        expected = study["outer_folds"] * len(config["quantiles"])
        progress.append({"study_id": study["study_id"], "completed": len(files), "expected": expected,
                         "percent": 100 * len(files) / expected, "status": state["status"]})
    report = Path(config["report_directory"])
    report.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(progress).to_csv(report / "progress.csv", index=False)
    if not pending and not running:
        break
    time.sleep(5)
(root / "status.json").write_text(json.dumps({"status": "completed", "fits": sum(row["expected"] for row in progress),
    "formal_validation_used": False, "manuscript_adoption": False}, indent=2) + "\n")
