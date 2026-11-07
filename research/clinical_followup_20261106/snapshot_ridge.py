import json
from pathlib import Path
import shutil
import sys

import pandas as pd
import yaml

from panel import study_config

study = study_config(sys.argv[1])
root = Path(study["output_root"])
reference = root/"implementation"
checks = pd.read_csv("results/preflight/clinical_followup_20261106/ridge_equivalence.csv")
assert checks.prediction_difference.max() < 1e-7
assert checks.lambda_difference.eq(0).all()
assert checks.solver.nunique() == 2
snapshot = root/"implementation_ridge_20261107"
sources = {"config.yml": reference/"config.yml", "run_discovery.R": reference/"run_discovery.R",
           "panel.R": reference/"panel.R", "ridge_direction.R": reference/"ridge_direction.R",
           "prediction_core.R": Path("research/empirical_selection_20261101/prediction_core.R")}
assert yaml.safe_load(sources["config.yml"].read_text()) == study
snapshot.mkdir(exist_ok=True)
for name, path in sources.items():
    destination = snapshot/name
    if destination.exists():
        assert destination.read_bytes() == path.read_bytes()
    else:
        shutil.copy2(path,destination)
receipt = {"purpose":"equivalent ridge solution in the smaller Gram matrix",
           "reference_implementation":str(reference), "current_implementation":str(snapshot),
           "maximum_verified_prediction_difference":float(checks.prediction_difference.max()),
           "regularization_grid_verified":True, "formal_validation_used":False}
(snapshot/"provenance.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(snapshot)
