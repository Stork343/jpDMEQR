import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

config = yaml.safe_load(Path("config/application_ra_westchina_20261105.yml").read_text())
source = Path(config["source_directory"])
output = Path(config["input_directory"])
assert not output.exists()
clinical = pd.read_excel(source / "41467_2025_62032_MOESM3_ESM.xlsx", header=2)
protein = pd.read_excel(source / "41467_2025_62032_MOESM4_ESM.xlsx", header=2)
assert clinical["Sample ID"].is_unique and protein.protein.is_unique
assert protein.columns.is_unique
ra = clinical.loc[clinical.Group.eq("RA")].copy()
assert ra["Patient ID"].notna().all()
assert ra.groupby("Patient ID").Sex.nunique().eq(1).all()
time_columns = [f"Time point {i}" for i in [1, 2, 3]]
assert ra[time_columns].notna().sum(axis=1).eq(1).all()
ra["visit"] = ra[time_columns].notna().to_numpy().argmax(axis=1) + 1
assert not ra[["Patient ID", "visit"]].duplicated().any()
assert set(ra["Sample ID"]).issubset(protein.columns)
rows = ra.loc[ra[config["response"]].notna()].copy()
rows = rows.rename(columns={"Sample ID": "sample_id", "Patient ID": "strain_id", config["response"]: "y"})
rows = rows.sort_values(["strain_id", "visit", "sample_id"])
assert rows.y.gt(0).all() and np.isfinite(rows.y).all()
participants = rows.groupby("strain_id", sort=True).first()[["Sex"]].reset_index()
participants["role"] = "validation"
participants["outer_fold"] = 0
rng = np.random.default_rng(config["role_seed"])
for sex, indices in participants.groupby("Sex", sort=True).groups.items():
    ordered = rng.permutation(np.asarray(indices))
    selected = ordered[:round(config["training_fraction"] * len(ordered))]
    participants.loc[selected, "role"] = "discovery"
    participants.loc[selected, "outer_fold"] = np.arange(len(selected)) % config["outer_folds"] + 1
rows = rows.merge(participants[["strain_id", "role", "outer_fold"]], on="strain_id", validate="many_to_one")
rows["condition"] = "das28"
ccp = rows["CCP(U/ml)"].astype("string").str.strip()
upper = ccp.eq(">500").fillna(False)
lower = ccp.eq("<17").fillna(False)
numeric_ccp = pd.to_numeric(ccp.mask(upper | lower), errors="raise")
acpa = numeric_ccp.ge(17).where(numeric_ccp.notna(), np.nan).astype(float)
acpa.loc[upper] = 1.0
acpa.loc[lower] = 0.0
clinical_values = pd.DataFrame({"sample_id": rows.sample_id, "C_age": rows["Age(Years)"],
    "C_male": rows.Sex.eq("Male").astype(float), "C_acpa": acpa,
    "C_visit2": rows.visit.eq(2).astype(float), "C_visit3": rows.visit.eq(3).astype(float)})
excluded = set(config["excluded_genes"])
protein["excluded_response_component"] = protein.gene.fillna("").map(lambda value: bool(set(value.split(";")) & excluded))
catalog = protein[["protein", "gene", "excluded_response_component"]].copy()
selected = protein.loc[~protein.excluded_response_component].set_index("protein")
P = selected[rows.sample_id].T
assert np.isfinite(P.to_numpy()[P.notna().to_numpy()]).all()
assert (P.to_numpy()[P.notna().to_numpy()] > 0).all()
P = np.log2(P)
P.columns = ["P_" + x for x in P.columns]
P.index.name = "sample_id"
features = clinical_values.set_index("sample_id").join(P, validate="one_to_one")
assert features.index.equals(pd.Index(rows.sample_id))
output.mkdir(parents=True)
participants.to_csv(output / "roles.csv", index=False)
catalog.to_csv(output / "protein_catalog.csv", index=False)
features.to_csv(output / "features.csv.gz")
columns = ["sample_id", "strain_id", "visit", "y", "condition", "role", "outer_fold"]
for role in ["discovery", "validation"]:
    rows.loc[rows.role.eq(role), columns].to_csv(output / f"{role}.csv", index=False)
report = Path(config["report_directory"])
report.mkdir(parents=True, exist_ok=True)
counts = rows.groupby("role").agg(patients=("strain_id", "nunique"), records=("sample_id", "size"))
counts.to_csv(report / "cohort_counts.csv")
cluster_sizes = rows.groupby("strain_id").size()
status = {"status": "prepared", "source_doi": config["source_doi"], "correction_doi": config["correction_doi"],
    "ra_source_patients": ra["Patient ID"].nunique(), "ra_source_records": len(ra),
    "patients": rows.strain_id.nunique(), "records": len(rows),
    "patients_with_repeats": int(cluster_sizes.ge(2).sum()), "maximum_cluster_size": int(cluster_sizes.max()),
    "protein_candidates": len(selected), "excluded_protein_components": int(catalog.excluded_response_component.sum()),
    "clinical_predictors": list(clinical_values.columns[1:]), "outcome_imputed": False,
    "formal_validation_used": False, "manuscript_adoption": False}
(report / "cohort_status.json").write_text(json.dumps(status, indent=2) + "\n")
(report / "data_manifest.json").write_text((source / "manifest.json").read_text())
print(json.dumps(status, indent=2))
print(counts.to_string())
