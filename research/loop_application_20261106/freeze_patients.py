import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd
import yaml

config = yaml.safe_load(Path("config/application_loop_20261106.yml").read_text())
source = Path(config["source_directory"])
inputs = Path(config["input_directory"])
report = Path(config["report_directory"])
inputs.mkdir(parents=True, exist_ok=True)
report.mkdir(parents=True, exist_ok=True)
if (source / "Loop.zip").exists():
    with ZipFile(source / "Loop.zip") as archive:
        csv = archive.read("Data Tables/PtRoster.txt")
    if (source / "PtRoster.txt").exists():
        assert (source / "PtRoster.txt").read_bytes() == csv
    else:
        (source / "PtRoster.txt").write_bytes(csv)
roster = pd.read_csv(source / "PtRoster.txt", sep="|", dtype={"PtID": str})
cohort = roster.loc[roster.PtCohort.eq("A"), ["PtID", "AgeAtEnrollment", "LoopUseTimeAtEnroll", "EnrollDt",
                                              "PtTimezoneOffset"]].copy()
assert cohort.PtID.is_unique and cohort[["PtID", "AgeAtEnrollment", "EnrollDt"]].notna().all().all()
cohort["pediatric"] = cohort.AgeAtEnrollment < 18
cohort["use_at_enrollment"] = cohort.LoopUseTimeAtEnroll.fillna("Unknown")
random = np.random.default_rng(config["role_seed"])
cohort = cohort.sort_values("PtID").set_index("PtID")
cohort["role"] = "discovery"
cohort["outer_fold"] = 0
for _, cell in cohort.groupby(["pediatric", "use_at_enrollment"], sort=True):
    ids = random.permutation(cell.index)
    reserved = int(np.rint(config["validation_fraction"] * len(ids)))
    cohort.loc[ids[:reserved], "role"] = "reserved"
    remaining = ids[reserved:]
    shift = int(random.integers(config["outer_folds"]))
    cohort.loc[remaining, "outer_fold"] = (np.arange(len(remaining)) + shift) % config["outer_folds"] + 1
assert cohort.loc[cohort.role.eq("reserved")].shape[0] >= config["minimum_reserved_roster"]
assert set(cohort.loc[cohort.role.eq("discovery"), "outer_fold"]) == set(range(1, config["outer_folds"] + 1))
csv = cohort.to_csv()
if (inputs / "roles.csv").exists():
    assert (inputs / "roles.csv").read_text() == csv
else:
    (inputs / "roles.csv").write_text(csv)
inventory = cohort.groupby(["role", "outer_fold", "pediatric", "use_at_enrollment"]).size()
inventory.rename("patients").to_csv(report / "role_inventory.csv")
status = {"source": "Jaeb Loop PtRoster cohort A", "registered": len(cohort),
          "discovery_patients": int(cohort.role.eq("discovery").sum()),
          "reserved_patients": int(cohort.role.eq("reserved").sum()),
          "formal_validation_used": False, "manuscript_adoption": False}
(report / "roles.json").write_text(json.dumps(status, indent=2) + "\n")
print(json.dumps(status))
