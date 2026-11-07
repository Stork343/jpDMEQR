import json
from pathlib import Path
import shutil
import sys

import yaml

from panel import study_config

study = study_config(sys.argv[1])
root = Path(study["output_root"])
reference = root/"implementation"
snapshot = root/"implementation_joint_20261107"
sources = {"config.yml":reference/"config.yml","panel.R":reference/"panel.R",
           "ridge_direction.R":reference/"ridge_direction.R",
           "run_discovery.R":Path("research/clinical_search_20261105/run_discovery.R"),
           "prediction_core.R":Path("research/empirical_selection_20261101/prediction_core.R"),
           "conic_search.R":Path("research/joint_profile_manuscript_20261026/conic_search.R")}
assert yaml.safe_load(sources["config.yml"].read_text()) == study
snapshot.mkdir(exist_ok=True)
for name,path in sources.items():
    target = snapshot/name
    if target.exists():
        assert target.read_bytes() == path.read_bytes()
    else:
        shutil.copy2(path,target)
state = {"solver":"clarabel_selected_qp","objective":"original joint profile with prediction ridge",
         "candidate_limit":"minimum of registered size and independent directions in all inner and final training matrices",
         "rank_tolerance":1e-5,"quantile_loss_and_working_penalties":study["working_penalties"],
         "formal_validation_used":False,"manuscript_adoption":False}
(snapshot/"provenance.json").write_text(json.dumps(state,indent=2)+"\n")
print(snapshot)
