import argparse
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import yaml

parser = argparse.ArgumentParser()
parser.add_argument("config", type=Path)
args = parser.parse_args()
config = yaml.safe_load(args.config.read_text())
root = Path(config["output_root"])
backup = Path("/Users/houjian/jpDMEQR-backups/empirical_selection_20261101/results/application") / root.name
backup.mkdir(parents=True, exist_ok=True)
progress = pd.read_csv(root / "progress.csv")
assert len(progress) == config["outer_folds"] * 15
completed = progress.loc[progress.status.eq("completed")]
for directory in completed.directory:
    source = Path(directory)
    status = json.loads((source / "status.json").read_text())
    assert status["status"] == "completed" and (source / "fit.rds").exists()
    destination = backup / source.relative_to(root).parent
    destination.mkdir(parents=True, exist_ok=True)
    subprocess.run(["rsync", "-a", str(source), str(destination) + "/"], check=True)
    for path in source.rglob("*"):
        if path.is_file():
            copied = destination / source.name / path.relative_to(source)
            assert copied.stat().st_size == path.stat().st_size
for source in [root / "implementation", *root.glob("outer_*/recipes.rds"),
               *root.glob("outer_*/training.csv"), *root.glob("outer_*/evaluation.csv"),
               root / "jobs.csv", root / "preparation.json", root / "preparation_summary.csv"]:
    assert source.exists()
    destination = backup / source.relative_to(root).parent
    destination.mkdir(parents=True, exist_ok=True)
    subprocess.run(["rsync", "-a", str(source), str(destination) + "/"], check=True)
status = {"study_id": config["study_id"], "completed_model_directories": len(completed),
          "registered_jobs": len(progress), "backup_directory": str(backup),
          "checked_utc": datetime.now(timezone.utc).isoformat(),
          "full_screen_backed_up": len(completed) == len(progress)}
(backup / "backup_status.json").write_text(json.dumps(status, indent=2) + "\n")
(Path(config["report_directory"]) / "backup_status.json").write_text(json.dumps(status, indent=2) + "\n")
print(json.dumps(status, indent=2))
