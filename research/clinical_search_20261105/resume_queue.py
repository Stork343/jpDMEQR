import argparse
import json
import os
import subprocess
import time
from pathlib import Path

import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("root", type=Path)
parser.add_argument("workers", type=int, choices=[1, 2, 3, 4])
args = parser.parse_args()
progress = pd.read_csv(args.root / "progress.csv")
pending = progress.loc[progress.status.eq("running")].copy()
while len(pending):
    completed = []
    for index, job in pending.iterrows():
        status_file = Path(job.directory) / "status.json"
        if status_file.exists():
            status = json.loads(status_file.read_text())
            assert status["status"] == "completed"
            completed.append(index)
        else:
            os.kill(int(job.pid), 0)
    pending = pending.drop(index=completed)
    if len(pending):
        time.sleep(2)
print(f"Continuing queue with {args.workers} workers", flush=True)
subprocess.run(["Rscript", "research/yeast_application_20261101/run_queue.R",
                str(args.root / "implementation"), str(args.workers), str(args.root)], check=True)
