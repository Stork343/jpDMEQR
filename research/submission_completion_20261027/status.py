import csv
import json
from collections import Counter
from pathlib import Path


def read_rows(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


reports = []
for root in sorted(Path("results/raw/simulation").glob("submission_2026*_main_*")):
    scheduled = read_rows(root / "scheduled_replications.csv")
    completed = sum(
        int(row["completed"])
        for path in root.glob("batch_*_progress.csv")
        for row in read_rows(path)
    )
    reports.append({
        "queue": root.name,
        "completed_datasets": completed,
        "registered_datasets": len(scheduled),
        "finished": (root / "completed.rds").exists(),
    })
for name in ("submission_20261027_v2", "submission_20261028_qgee"):
    root = Path("results/raw/application_comparison") / name
    if not (root / "status.csv").exists():
        continue
    rows = read_rows(root / "status.csv")
    states = Counter(row["state"] for row in rows)
    reports.append({
        "queue": name,
        "successful_jobs": sum(row["state"] == "completed" and row["exit_status"] == "0" for row in rows),
        "failed_jobs": sum(row["state"] == "completed" and row["exit_status"] != "0" for row in rows),
        "registered_jobs": len(rows),
        "states": dict(states),
        "finished": (root / "completed.rds").exists(),
    })
print(json.dumps(reports, indent=2))
