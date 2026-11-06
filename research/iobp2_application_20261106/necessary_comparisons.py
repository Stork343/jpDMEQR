from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import os
import subprocess
import sys

import numpy as np
import pandas as pd
import yaml

config = yaml.safe_load(Path(sys.argv[1]).read_text())
root = Path(config["output_root"])
report = Path(config["report_directory"])
assert json.loads((report / "eligibility.json").read_text())["passed"]
implementation = root / "implementation"
assert yaml.safe_load((implementation / "config.yml").read_text()) == config
jobs = pd.read_csv(root / "jobs.csv").query("method in ['CLINICAL_QR', 'ENET']")
assert len(jobs) == config["outer_folds"] * 2
environment = dict(os.environ, TMPDIR=str(Path("results/scratch/runtime").resolve()), OMP_NUM_THREADS="1",
                   OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1", R_MAX_VSIZE="24Gb")


def fit(job):
    directory = Path(job.directory)
    exit_code = None
    if not (directory / "status.json").exists():
        with (root / f"{job.trait}_{job.method}.log").open("a") as stream:
            process = subprocess.run(["Rscript", str(implementation / "run_discovery.R"), job.trait, job.method,
                            job.run_id, "NA", str(implementation.resolve())], env=environment,
                           stdout=stream, stderr=subprocess.STDOUT, check=True)
            exit_code = process.returncode
    status = json.loads((directory / "status.json").read_text())
    assert status["status"] == "completed" and not status["formal_validation_used"]
    predictions = pd.read_csv(directory / "prediction.csv")
    assert predictions.method.eq(job.method).all()
    rows = pd.read_csv(root / job.trait / "evaluation.csv")
    assert rows.role.eq("discovery_evaluation").all() and len(predictions) == len(rows) * len(config["quantiles"])
    frame = predictions.merge(rows[["sample_id", "strain_id", "y", "outer_fold"]], on=["sample_id", "strain_id"],
                              validate="many_to_one")
    assert len(frame) == len(predictions) and set(frame.sample_id) == set(rows.sample_id)
    print(f"{job.trait} {job.method} completed", flush=True)
    return frame, {"trait": job.trait, "method": job.method, "directory": str(directory),
                   "status": status["status"], "exit_code": exit_code}


with ThreadPoolExecutor(max_workers=2) as executor:
    results = list(executor.map(fit, jobs.itertuples()))
records = pd.concat([frame for frame, _ in results], ignore_index=True)
progress = pd.DataFrame([state for _, state in results])
progress_file = root / "necessary_progress.csv"
if progress_file.exists():
    previous = pd.read_csv(progress_file).set_index(["trait", "method"]).exit_code
    progress["exit_code"] = progress.set_index(["trait", "method"]).exit_code.fillna(previous).to_numpy()
progress.to_csv(progress_file, index=False)
groups = ["method", "sample_id", "strain_id", "y", "outer_fold"]
wide = records.pivot(index=groups, columns="tau", values="prediction")
assert list(wide.columns) == config["quantiles"] and np.isfinite(wide).all().all()
sorted_predictions = wide.copy()
sorted_predictions.iloc[:, :] = np.maximum(0, np.sort(wide.to_numpy(), axis=1))
ordered = sorted_predictions.stack().rename("prediction").reset_index()
results = []
for scoring, frame in [("raw", records), ("primary", ordered)]:
    residual = frame.y - frame.prediction
    loss = np.maximum(frame.tau * residual, (frame.tau - 1) * residual)
    patient = frame.assign(loss=loss).groupby(["method", "tau", "strain_id"]).loss.mean()
    summary = patient.groupby(["method", "tau"]).agg(patients="size", loss="mean", sd="std").reset_index()
    summary["se"] = summary.sd / np.sqrt(summary.patients)
    summary["scoring"] = scoring
    results.append(summary)
pd.concat(results).to_csv(report / "discovery_loss.csv", index=False)
