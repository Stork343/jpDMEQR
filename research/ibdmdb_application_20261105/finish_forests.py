import json
import subprocess
import time
from pathlib import Path

import pandas as pd

root = Path("results/raw/application/ibdmdb_forests_20261105")
while True:
    jobs = pd.read_csv(root / "progress.csv")
    assert len(jobs) == 10
    if jobs.status.eq("failed").any():
        raise RuntimeError("计算清单中存在失败记录")
    if jobs.status.eq("completed").all():
        break
    time.sleep(15)
with Path("results/logs/yeast_1011_20261104/ibdmdb_forests_report.log").open("w") as log:
    subprocess.run(["Rscript", "research/ibdmdb_application_20261105/report.R",
                    "results/raw/application/ibdmdb_median_20261105",
                    "results/preflight/empirical_selection_20261101/ibdmdb_forests", "0.5", str(root)],
                   stdout=log, stderr=subprocess.STDOUT, check=True)
subprocess.run(["rsync", "-a", str(root),
                "/Users/houjian/jpDMEQR-backups/empirical_selection_20261101/results/application/"], check=True)
print(json.dumps({"status": "completed", "formal_validation_used": False, "manuscript_adoption": False}), flush=True)
