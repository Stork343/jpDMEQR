import json
import subprocess
import time
from pathlib import Path

import pandas as pd

queue_root = Path("results/raw/application/yeast_stress_20261104/screening_20261104")
progress = queue_root / "progress.csv"
logs = Path("results/logs/yeast_1011_20261104")
while True:
    jobs = pd.read_csv(progress)
    assert len(jobs) == 12
    if jobs.status.eq("failed").any():
        raise RuntimeError("计算清单中存在失败记录")
    if jobs.status.eq("completed").all():
        break
    time.sleep(15)
with (logs / "stress_screening_report.log").open("w") as log:
    subprocess.run(["Rscript", "research/yeast_stress_application_20261104/report_screening.R"],
                   stdout=log, stderr=subprocess.STDOUT, check=True)
comparison = pd.read_csv("results/preflight/empirical_selection_20261101/yeast_stress_prediction/comparisons.csv")
eligible = comparison.loc[comparison.exact_comparison_required, "assay"].tolist()
directory = Path("results/raw/application/yeast_stress_20261104/exact_comparisons_20261104")
assert not directory.exists()
directory.mkdir()
if eligible:
    followup = pd.DataFrame({"trait": eligible, "method": "POOLED", "tau": 0.25, "run_id": "tau025"})
    followup["job_id"] = range(1, len(followup) + 1)
    followup["directory"] = "results/raw/application/yeast_stress_20261104/discovery/" + followup.trait + "/POOLED_tau025"
    followup["status"] = "pending"
    followup.to_csv(directory / "jobs.csv", index=False)
    with (logs / "stress_exact_queue.log").open("w") as log:
        subprocess.run(["Rscript", "research/yeast_application_20261101/run_queue.R",
                        str(queue_root / "implementation"), "2", str(directory)],
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    with (logs / "stress_exact_report.log").open("w") as log:
        subprocess.run(["Rscript", "research/yeast_stress_application_20261104/report_exact.R"],
                       stdout=log, stderr=subprocess.STDOUT, check=True)
subprocess.run([
    "rsync", "-a", "results/raw/application/yeast_stress_20261104",
    "/Users/houjian/jpDMEQR-backups/empirical_selection_20261101/results/application/",
], check=True)
(directory / "stage_complete.json").write_text(json.dumps({
    "status": "completed", "exact_comparison_assays": eligible,
    "formal_validation_used": False, "manuscript_adoption": False,
}, indent=2))
print(json.dumps({"screening": "completed", "exact_comparisons": eligible}), flush=True)
