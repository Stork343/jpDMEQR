import argparse
import json
import os
import time
from pathlib import Path

import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("root", type=Path)
args = parser.parse_args()
progress = pd.read_csv(args.root / "progress.csv")
active = progress.loc[progress.status.eq("running")].copy()
while len(active):
    completed = []
    for index, job in active.iterrows():
        status_file = Path(job.directory) / "status.json"
        if status_file.exists():
            status = json.loads(status_file.read_text())
            assert status["status"] == "completed"
            assert (Path(job.directory) / "fit.rds").exists()
            assert (Path(job.directory) / "prediction.csv").exists()
            progress.loc[index, "status"] = "completed"
            completed.append(index)
        else:
            os.kill(int(job.pid), 0)
    active = active.drop(index=completed)
    if len(active):
        time.sleep(2)
progress.loc[progress.status.eq("pending"), "status"] = "not_run"
progress.to_csv(args.root / "progress.csv", index=False)
status = {"status": "closed", "analysis_role": "adjacent_visit_prediction_diagnostic",
          "registered_jobs": len(progress), "completed_jobs": int(progress.status.eq("completed").sum()),
          "unrun_jobs": int(progress.status.eq("not_run").sum()), "formal_validation_used": False,
          "manuscript_adoption": False}
(args.root / "execution_status.json").write_text(json.dumps(status, indent=2) + "\n")
print(json.dumps(status, indent=2))
