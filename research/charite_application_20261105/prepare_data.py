import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

config = yaml.safe_load(Path("config/application_charite_hb_20261105.yml").read_text())
source = Path(config["source_directory"])
output = Path(config["input_directory"])
assert not output.exists()
for filename in ["mmc7.xlsx", "mmc8.xlsx", "mmc9.xlsx"]:
    with zipfile.ZipFile(source / filename) as archive:
        assert archive.testzip() is None
values = pd.read_excel(source / "mmc7.xlsx").set_index("Measurement")
metadata = pd.read_excel(source / "mmc8.xlsx").set_index("File")
annotation = pd.read_excel(source / "mmc9.xlsx")
assert values.index.is_unique and values.columns.is_unique and metadata.index.is_unique
assert set(values.columns) == set(metadata.index)
assert annotation.Genes.is_unique
metadata = metadata.loc[values.columns].copy()
metadata["sample_id"] = metadata.index
metadata["strain_id"] = metadata["Patient Id"].map(lambda x: f"P{int(x):04d}")
metadata["y"] = pd.to_numeric(values.loc[config["response"]], errors="raise").to_numpy()
cohort = metadata.loc[metadata.y.notna()].copy()
assert cohort.y.gt(0).all() and np.isfinite(cohort.y).all()
assert not cohort[["strain_id", "Time point - days since admission"]].duplicated().any()
cohort = cohort.rename(columns={"Age": "age", "BMI": "bmi", "Charlson index": "charlson",
                               "WHO grade": "who", "Time point - days since admission": "admission_day"})
assert cohort[["age", "charlson", "who", "admission_day"]].notna().all().all()
assert cohort.admission_day.ge(0).all() and cohort.who.isin([3, 4, 5, 6, 7]).all()
cohort = cohort.sort_values(["strain_id", "admission_day", "sample_id"])
participants = cohort.groupby("strain_id", sort=True).first()[["who", "age"]].reset_index()
participants["stratum"] = np.where(participants.who.ge(6), "invasive", "noninvasive")
participants["role"] = "validation"
participants["inner_fold"] = 0
rng = np.random.default_rng(config["role_seed"])
for stratum, indices in participants.groupby("stratum", sort=True).groups.items():
    order = rng.permutation(np.asarray(indices))
    ntrain = int(round(config["training_fraction"] * len(order)))
    selected = order[:ntrain]
    participants.loc[selected, "role"] = "training"
    participants.loc[selected, "inner_fold"] = np.arange(ntrain) % config["inner_folds"] + 1
cohort = cohort.merge(participants[["strain_id", "role", "inner_fold"]], on="strain_id", validate="many_to_one")
cohort["condition"] = "hemoglobin"
keep = ["sample_id", "strain_id", "admission_day", "age", "bmi", "charlson", "who", "y", "condition", "role", "inner_fold"]
cohort = cohort[keep].sort_values(["strain_id", "admission_day", "sample_id"])
excluded = set(config["excluded_gene_tokens"])
catalog = annotation.loc[annotation.Genes.isin(values.index) & annotation.Peptides.ge(config["protein_minimum_peptides"])].copy()
catalog["excluded_response_protein"] = catalog.Genes.map(lambda genes: bool(set(genes.split(";")) & excluded))
proteins = catalog.loc[~catalog.excluded_response_protein, "Genes"]
P = values.loc[proteins, cohort.sample_id].T
assert np.isfinite(P.to_numpy()[P.notna().to_numpy()]).all()
assert (P.to_numpy()[P.notna().to_numpy()] > 0).all()
output.mkdir(parents=True)
participants.to_csv(output / "roles.csv", index=False)
catalog.to_csv(output / "protein_catalog.csv", index=False)
P.to_csv(output / "protein_abundance.csv.gz", index_label="sample_id")
for role in ["training", "validation"]:
    cohort.loc[cohort.role.eq(role)].to_csv(output / f"{role}.csv", index=False)
report = Path(config["report_directory"])
report.mkdir(parents=True, exist_ok=True)
counts = cohort.groupby("role").agg(patients=("strain_id", "nunique"), records=("sample_id", "size"),
                                    bmi_missing=("bmi", lambda x: int(x.isna().sum())))
counts.to_csv(report / "cohort_counts.csv")
summary = {"status": "prepared", "source_doi": config["source_doi"], "source_rows": len(metadata),
           "patients": len(participants), "records": len(cohort), "protein_candidates": len(proteins),
           "excluded_response_protein_groups": int(catalog.excluded_response_protein.sum()),
           "maximum_cluster_size": int(cohort.groupby("strain_id").size().max()),
           "protein_missing_fraction": float(P.isna().to_numpy().mean()),
           "outcome_imputed": False, "response_unit": config["response_unit"],
           "formal_validation_used": False, "outcome_associations_fitted": False,
           "manuscript_adoption": False}
(report / "cohort_status.json").write_text(json.dumps(summary, indent=2) + "\n")
print(counts.to_string())
print(json.dumps(summary, indent=2))
