from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import subprocess
import sys

import yaml

python = sys.executable
source = Path("research/loop_application_20261106")
common = Path("research/iobp2_application_20261106")
smooth = Path("research/clinical_quantile_20261105")
clinical = Path("research/clinical_search_20261105")
study_path = Path("config/application_loop_20261106.yml")
method_path = Path("config/clinical_quantile_loop_20261106.yml")
study = yaml.safe_load(study_path.read_text())
root = Path(study["output_root"])
report = Path(study["report_directory"])
inputs = Path(study["input_directory"])
environment = dict(os.environ, TMPDIR=str(Path("results/scratch/runtime").resolve()), OMP_NUM_THREADS="1",
                   OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1", R_MAX_VSIZE="24Gb",
                   PYTHONDONTWRITEBYTECODE="1")


def run(command):
    print(" ".join(map(str, command)), flush=True)
    subprocess.run(list(map(str, command)), check=True, env=environment)


if not (Path(study["source_directory"]) / "verified.json").exists():
    run([python, source / "download.py"])
run([python, source / "freeze_patients.py"])
if not (report / "source_counts.json").exists():
    run([python, source / "extract_cgm.py"])
run([python, common / "build_discovery.py", study_path])
run([python, common / "check_panel.py", study_path])
if not json.loads((report / "eligibility.json").read_text())["passed"]:
    raise SystemExit(0)
if not (root / "preparation.json").exists():
    run(["Rscript", clinical / "prepare_discovery.R", study_path])
with ThreadPoolExecutor(max_workers=2) as executor:
    baseline = executor.submit(run, [python, common / "necessary_comparisons.py", study_path])
    model = executor.submit(run, [python, smooth / "run_smooth_queue.py", method_path])
    baseline.result()
    model.result()
run(["Rscript", clinical / "verify_screen.R", study_path, "necessary_progress.csv"])
run([python, smooth / "summarise.py", method_path])
run(["Rscript", smooth / "verify_predictions.R", method_path])
run([python, common / "report_necessary.py", method_path])
