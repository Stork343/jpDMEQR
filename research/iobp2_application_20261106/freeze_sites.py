from itertools import combinations
from pathlib import Path
import json
from zipfile import ZipFile

import numpy as np
import pandas as pd
import yaml

config = yaml.safe_load(Path("config/application_iobp2_20261106.yml").read_text())
destination = Path(config["input_directory"])
report = Path(config["report_directory"])
destination.mkdir(parents=True, exist_ok=True)
report.mkdir(parents=True, exist_ok=True)
with ZipFile(Path(config["source_directory"]) / "IOBP2_RCT.zip") as archive:
    with archive.open("Data Tables/IOBP2PtRoster.txt") as stream:
        roster = pd.read_csv(stream, sep="|", dtype={"PtID": str})
roster = roster.loc[roster.TrtGroup.isin(["BP", "Control"]),
                    ["PtID", "SiteID", "TrtGroup", "AgeAsofEnrollDt", "RandDt"]].copy()
assert len(roster) == 326 and roster.PtID.is_unique and roster.notna().all().all()
roster["pediatric"] = roster.AgeAsofEnrollDt < 18
strata = pd.crosstab(roster.SiteID, [roster.TrtGroup, roster.pediatric])
sites = np.random.default_rng(config["role_seed"]).permutation(sorted(strata.index))
counts = strata.loc[sites].to_numpy()
totals = counts.sum(axis=0)
fraction = config["validation_fraction"]
best_key = None
for number in range(4, len(sites) - 6 + 1):
    subsets = np.asarray(list(combinations(range(len(sites)), number)))
    distances = ((counts[subsets].sum(axis=1) - fraction * totals) ** 2 / totals).sum(axis=1)
    index = int(distances.argmin())
    key = (float(distances[index]), abs(number - fraction * len(sites)), tuple(subsets[index]))
    if best_key is None or key < best_key:
        best_key = key
        reserved = sites[subsets[index]]

fold_sizes = np.zeros(config["outer_folds"], dtype=int)
folds = {}
site_sizes = strata.sum(axis=1)
for site in sorted((site for site in sites if site not in reserved), key=lambda site: -site_sizes[site]):
    fold = int(fold_sizes.argmin())
    folds[site] = fold + 1
    fold_sizes[fold] += site_sizes[site]
roster["role"] = np.where(roster.SiteID.isin(reserved), "reserved", "discovery")
roster["outer_fold"] = roster.SiteID.map(folds).fillna(0).astype(int)
assert roster.groupby("SiteID").role.nunique().eq(1).all()
assert roster.groupby("SiteID").outer_fold.nunique().eq(1).all()
assert roster.role.eq("reserved").sum() >= config["minimum_reserved_roster"]
roles_path = destination / "roles.csv"
csv = roster.sort_values("PtID").to_csv(index=False)
if roles_path.exists():
    assert roles_path.read_text() == csv
else:
    roles_path.write_text(csv)
inventory = roster.groupby(["role", "SiteID", "outer_fold", "TrtGroup", "pediatric"]).size()
inventory.rename("randomized_patients").to_csv(report / "role_inventory.csv")
metadata = {"roster": len(roster), "archive_site_codes": len(sites),
            "discovery_patients": int(fold_sizes.sum()), "reserved_patients": int(roster.role.eq("reserved").sum()),
            "discovery_site_codes": len(folds), "reserved_site_codes": len(reserved),
            "discovery_patients_per_fold": fold_sizes.tolist(), "source": "IOBP2PtRoster",
            "formal_validation_used": False, "manuscript_adoption": False}
(report / "roles.json").write_text(json.dumps(metadata, indent=2) + "\n")
print(json.dumps(metadata))
