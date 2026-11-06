import json
import os
from pathlib import Path
import subprocess
import sys

import pandas as pd
import yaml

from panel import settings_directory

name = sys.argv[1]
settings = settings_directory / f"{name}.yml"
method_path = settings_directory / f"cq_{name}.yml"
config = yaml.safe_load(settings.read_text())
method = yaml.safe_load(method_path.read_text())
root = Path(config["output_root"])
report = Path(config["report_directory"])
previous = json.loads((report / "decision.json").read_text())
assert previous["necessary_pass"] and not previous["formal_validation_used"]
if config.get("maximum_high_features"):
    assert json.loads((report / "gene_selection_verification.json").read_text())["status"] == "passed"
environment = dict(os.environ, TMPDIR=str(Path("results/scratch/runtime").resolve()), OMP_NUM_THREADS="1",
    OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1", PYTHONDONTWRITEBYTECODE="1", R_MAX_VSIZE="24Gb")
commands = [
    ["Rscript", "research/yeast_application_20261101/run_queue.R", root / "implementation", "2", root],
    [sys.executable, "research/clinical_search_20261105/summarise_discovery.py", settings, "--comparators-only"],
    ["Rscript", "research/clinical_search_20261105/verify_screen.R", settings],
    [sys.executable, "research/clinical_quantile_20261105/summarise.py", method_path],
    ["Rscript", "research/clinical_quantile_20261105/verify_predictions.R", method_path]]
for command in commands:
    subprocess.run(list(map(str, command)), check=True, env=environment)
comparison = pd.read_csv(Path(method["report_directory"]) / "discovery_comparisons.csv").query("decision_quantile")
assert set(comparison.comparator) == {"CLINICAL_QR", "ENET", "PRESS_RIDGE", "POOLED", "POOLED_OFFSET", "JOINT_DIRECT", "JOINT_OFFSET"}
gains = comparison.set_index("comparator").relative_improvement
passed = bool(gains.min() >= config["screening_minimum_relative_gain"])
state = {"study_id": config["study_id"], "comparators": gains.to_dict(), "full_linear_pass": passed,
         "phase": "nonlinear_comparisons" if passed else "closed_after_full_linear_comparisons",
         "formal_validation_used": False, "manuscript_adoption": False}
(report / "full_linear_decision.json").write_text(json.dumps(state, indent=2) + "\n")
print(json.dumps(state))
