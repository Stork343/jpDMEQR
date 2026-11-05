import json
from pathlib import Path

import pandas as pd
import yaml

config = yaml.safe_load(Path("config/clinical_quantile_20261105.yml").read_text())
output = Path("results/preflight/clinical_quantile_20261105")
rows = []
for path in config["source_studies"]:
    study = yaml.safe_load(Path(path).read_text())
    root = Path(config["output_root"]) / study["study_id"]
    completed = 0
    for fold in range(1, study["outer_folds"] + 1):
        for tau in config["quantiles"]:
            directory = root / f"outer_{fold}" / f"tau{round(100 * tau):03d}"
            status_path = directory / "status.json"
            if status_path.exists():
                state = json.loads(status_path.read_text())
                assert state["status"] == "completed"
                assert (directory / "fit.rds").is_file() and (directory / "prediction.csv").is_file()
                completed += 1
    expected = study["outer_folds"] * len(config["quantiles"])
    state = json.loads((root / "status.json").read_text()) if root.exists() else {"status": "pending"}
    rows.append({"study_id": study["study_id"], "completed": completed, "expected": expected,
                 "percent": 100 * completed / expected, "status": state["status"]})
table = pd.DataFrame(rows)
table.to_csv(output / "progress.csv", index=False)
print(table.to_string(index=False))
