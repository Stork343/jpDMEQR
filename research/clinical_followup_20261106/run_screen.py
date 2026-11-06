from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import subprocess
import sys

import pandas as pd
import yaml

from panel import settings_directory

name = sys.argv[1]
config = settings_directory / f"{name}.yml"
method = settings_directory / f"cq_{name}.yml"
study = yaml.safe_load(config.read_text())
report = Path(study["report_directory"])
root = Path(study["output_root"])
smooth = Path("research/clinical_quantile_20261105")
clinical = Path("research/clinical_search_20261105")
helper = Path("research/iobp2_application_20261106")
environment = dict(os.environ, TMPDIR=str(Path("results/scratch/runtime").resolve()), OMP_NUM_THREADS="1",
                   OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1", NUMBA_NUM_THREADS="1",
                   R_MAX_VSIZE="24Gb", PYTHONDONTWRITEBYTECODE="1")


def run(command):
    print(" ".join(map(str, command)), flush=True)
    subprocess.run(list(map(str, command)), check=True, env=environment)


eligibility = json.loads((report / "eligibility.json").read_text())
if not eligibility["passed"]:
    print(json.dumps(eligibility))
    raise SystemExit(0)
if not (root / "preparation.json").exists():
    run(["Rscript", clinical / "prepare_discovery.R", config])
with ThreadPoolExecutor(max_workers=2) as executor:
    baselines = executor.submit(run, [sys.executable, helper / "necessary_comparisons.py", config])
    candidate = executor.submit(run, [sys.executable, smooth / "run_smooth_queue.py", method])
    baselines.result()
    candidate.result()
run(["Rscript", clinical / "verify_screen.R", config, "necessary_progress.csv"])
run([sys.executable, smooth / "summarise.py", method])
run(["Rscript", smooth / "verify_predictions.R", method])
method_config = yaml.safe_load(method.read_text())
result = Path(method_config["report_directory"])
comparison = pd.read_csv(result / "discovery_comparisons.csv").query("decision_quantile")
assert set(comparison.comparator) == {"CLINICAL_QR", "ENET"}
assert json.loads((report / "screen_audit.json").read_text())["status"] == "passed"
assert json.loads((result / "prediction_acceptance.json").read_text())["status"] == "passed"
gain = comparison.set_index("comparator").relative_improvement
passed = bool(gain.min() >= study["screening_minimum_relative_gain"])
state = {"study_id": study["study_id"], "patients": eligibility["discovery_patients"], "records": eligibility["records"],
         "decision_quantile": study["decision_quantile"], "relative_gains": gain.to_dict(), "necessary_pass": passed,
         "phase": "additional_comparisons" if passed else "closed_after_necessary_comparisons",
         "formal_validation_used": False, "manuscript_adoption": False}
(report / "decision.json").write_text(json.dumps(state, indent=2) + "\n")
print(json.dumps(state))
