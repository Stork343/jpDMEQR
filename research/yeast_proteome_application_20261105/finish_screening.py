import json
import subprocess
import time
from pathlib import Path

import pandas as pd

queue_root = Path("results/raw/application/yeast_proteome_20261105/screening_20261105")
logs = Path("results/logs/yeast_1011_20261104")
while True:
    jobs = pd.read_csv(queue_root / "progress.csv")
    assert len(jobs) == 14
    if jobs.status.eq("failed").any():
        raise RuntimeError("计算清单中存在失败记录")
    if jobs.status.eq("completed").all():
        break
    time.sleep(15)
with (logs / "proteome_screening_report.log").open("w") as log:
    subprocess.run(["Rscript", "research/yeast_proteome_application_20261105/report.R", "screening"],
                   stdout=log, stderr=subprocess.STDOUT, check=True)
comparison = pd.read_csv("results/preflight/empirical_selection_20261101/yeast_proteome_screening/comparisons.csv")
eligible = comparison.loc[comparison.relative_gain.gt(0) & comparison.comparator.eq("ENET"), "assay"].tolist()
directory = Path("results/raw/application/yeast_proteome_20261105/exact_comparisons_20261105")
assert not directory.exists()
directory.mkdir()
if eligible:
    followup = pd.DataFrame({"trait": eligible, "method": "POOLED", "tau": 0.25, "run_id": "tau025"})
    followup["job_id"] = range(1, len(followup) + 1)
    followup["directory"] = "results/raw/application/yeast_proteome_20261105/discovery/" + followup.trait + "/POOLED_tau025"
    followup["status"] = "pending"
    followup.to_csv(directory / "jobs.csv", index=False)
    with (logs / "proteome_exact_queue.log").open("w") as log:
        subprocess.run(["Rscript", "research/yeast_application_20261101/run_queue.R",
                        str(queue_root / "implementation"), "2", str(directory)],
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    with (logs / "proteome_exact_report.log").open("w") as log:
        subprocess.run(["Rscript", "research/yeast_proteome_application_20261105/report.R", "exact"],
                       stdout=log, stderr=subprocess.STDOUT, check=True)
subprocess.run([
    "rsync", "-a", "results/raw/application/yeast_proteome_20261105",
    "/Users/houjian/jpDMEQR-backups/empirical_selection_20261101/results/application/",
], check=True)
(directory / "stage_complete.json").write_text(json.dumps({
    "status": "completed", "exact_comparison_assays": eligible,
    "formal_validation_used": False, "manuscript_adoption": False,
}, indent=2))
print(json.dumps({"screening": "completed", "exact_comparisons": eligible}), flush=True)
