import json
import subprocess
import time
from pathlib import Path

import pandas as pd

directory = Path("results/raw/application/arabidopsis_20261105/screening_20261105")
while True:
    jobs = pd.read_csv(directory / "progress.csv")
    assert len(jobs) == 5
    if jobs.status.eq("failed").any():
        raise RuntimeError("计算清单中存在失败记录")
    if jobs.status.eq("completed").all():
        break
    time.sleep(15)
with Path("results/logs/arabidopsis_20261105/comparison_report.log").open("w") as log:
    subprocess.run(["Rscript", "research/arabidopsis_application_20261105/report.R"],
                   stdout=log, stderr=subprocess.STDOUT, check=True)
subprocess.run([
    "rsync", "-a", "results/raw/application/arabidopsis_20261105",
    "/Users/houjian/jpDMEQR-backups/empirical_selection_20261101/results/application/",
], check=True)
print(json.dumps({"status": "completed", "formal_validation_used": False, "manuscript_adoption": False}), flush=True)
