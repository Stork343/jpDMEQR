import json
import subprocess
import time
from pathlib import Path

import pandas as pd

root = Path("results/raw/application/ridge_offset_forests_20261105")
output = Path("results/preflight/empirical_selection_20261101/ridge_offset_forests")
assays = ["Benomyl", "Caffeine"]
completed = set()
while len(completed) < len(assays):
    jobs = pd.read_csv(root / "progress.csv")
    assert len(jobs) == 4
    if jobs.status.eq("failed").any():
        raise RuntimeError("计算清单中存在失败记录")
    for assay in assays:
        if assay in completed:
            continue
        current = jobs.loc[jobs.trait.eq(assay)]
        assert len(current) == 2
        if not current.status.eq("completed").all():
            continue
        log_path = Path("results/logs/yeast_1011_20261104") / f"forest_{assay}_report.log"
        with log_path.open("w") as log:
            subprocess.run(["Rscript", "research/ridge_direction_20261105/report_forests.R", assay],
                           stdout=log, stderr=subprocess.STDOUT, check=True)
        destination = Path("/Users/houjian/jpDMEQR-backups/empirical_selection_20261101/results/application") / root.name
        destination.mkdir(parents=True, exist_ok=True)
        subprocess.run(["rsync", "-a", str(root / assay), str(destination) + "/"], check=True)
        completed.add(assay)
        print(json.dumps({"assay": assay, "status": "reported_and_backed_up"}), flush=True)
    if len(completed) < len(assays):
        time.sleep(15)

for filename in ["summary.csv", "comparisons.csv", "artifact_checks.csv", "cv_checks.csv",
                 "forest_settings.csv", "quantile_crossings.csv", "condition_metrics.csv", "screening_decision.csv"]:
    pd.concat([pd.read_csv(output / assay / filename) for assay in assays], ignore_index=True).to_csv(output / filename, index=False)
status = {"status": "completed", "phase": "method_development", "assays": assays,
          "formal_validation_used": False, "manuscript_adoption": False}
(output / "status.json").write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")
subprocess.run(["rsync", "-a", str(root),
                "/Users/houjian/jpDMEQR-backups/empirical_selection_20261101/results/application/"], check=True)
print(json.dumps(status), flush=True)
