import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

config = yaml.safe_load(Path("config/application_pregnancy_detroit_20261105.yml").read_text())
source = Path(config["source_directory"])
output = Path(config["input_directory"])
assert not output.exists()
path = source / "pone.0181468.s001.csv"
annotation = pd.read_csv(path, header=None, nrows=2, dtype=str)
frame = pd.read_csv(path, header=3)
assert list(frame.columns[:3]) == ["ID", "GA", "LatePE"] and frame.columns.is_unique
assert frame.groupby("ID").LatePE.nunique().eq(1).all()
assert not frame[["ID", "GA"]].duplicated().any()
assert frame.GA.gt(0).all() and np.isfinite(frame.GA).all()
frame["sample_id"] = [f"D{i:04d}" for i in range(1, len(frame) + 1)]
frame["strain_id"] = frame.ID.map(lambda value: f"P{int(value):03d}")
participants = frame.groupby("strain_id", sort=True).first()[["LatePE"]].reset_index()
participants["role"] = "validation"
participants["outer_fold"] = 0
rng = np.random.default_rng(config["role_seed"])
for group, indices in participants.groupby("LatePE", sort=True).groups.items():
    ordered = rng.permutation(np.asarray(indices))
    selected = ordered[:round(config["training_fraction"] * len(ordered))]
    participants.loc[selected, "role"] = "discovery"
    participants.loc[selected, "outer_fold"] = np.arange(len(selected)) % config["outer_folds"] + 1
rows = frame[["sample_id", "strain_id", "GA"]].rename(columns={"GA": "y"})
rows = rows.merge(participants[["strain_id", "role", "outer_fold"]], on="strain_id", validate="many_to_one")
rows["condition"] = "gestational_age"
proteins = [column for column in frame.columns if column.startswith("SL")]
assert len(proteins) == 1125
P = frame[proteins].to_numpy(dtype=float)
assert np.isfinite(P).all() and (P > 0).all()
features = pd.DataFrame(np.log2(P), columns=["P_" + column for column in proteins], index=frame.sample_id)
features.index.name = "sample_id"
catalog = pd.DataFrame({"protein": proteins, "symbol": annotation.iloc[0, 3:].to_numpy(),
                        "description": annotation.iloc[1, 3:].to_numpy()})
output.mkdir(parents=True)
participants.to_csv(output / "roles.csv", index=False)
catalog.to_csv(output / "protein_catalog.csv", index=False)
features.to_csv(output / "features.csv.gz")
for role in ["discovery", "validation"]:
    rows.loc[rows.role.eq(role)].to_csv(output / f"{role}.csv", index=False)
report = Path(config["report_directory"])
report.mkdir(parents=True, exist_ok=True)
counts = rows.groupby("role").agg(patients=("strain_id", "nunique"), records=("sample_id", "size"))
counts.to_csv(report / "cohort_counts.csv")
cluster_sizes = rows.groupby("strain_id").size()
status = {"status": "prepared", "source_doi": config["source_doi"],
    "patients": rows.strain_id.nunique(), "records": len(rows),
    "patients_with_repeats": int(cluster_sizes.ge(2).sum()), "maximum_cluster_size": int(cluster_sizes.max()),
    "protein_candidates": len(proteins), "clinical_predictors": [], "outcome_imputed": False,
    "source_protein_processing": "supplier_normalization_and_published_abundance_thresholding",
    "formal_validation_used": False, "manuscript_adoption": False}
(report / "cohort_status.json").write_text(json.dumps(status, indent=2) + "\n")
(report / "data_manifest.json").write_text((source / "manifest.json").read_text())
print(json.dumps(status, indent=2))
print(counts.to_string())
