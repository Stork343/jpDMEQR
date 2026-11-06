import csv
import gzip
from pathlib import Path

import numpy as np
import pandas as pd

from panel import finish_panel, study_config

study = study_config("ispy2_volume")
source = Path(study["source_directory"])
inputs = Path(study["input_directory"])
roles = pd.read_csv(inputs / "roles.csv", dtype={"patient": str}).set_index("patient")
cohort = roles.loc[roles.role.eq("discovery")].copy()
file = source / "expression.tsv.gz"
assert file.with_name(file.name + ".verified.json").exists()
with gzip.open(file, "rt") as stream:
    ids = next(csv.reader([stream.readline()], delimiter="\t"))
with gzip.open(source / "gene_header.part", "rt") as stream:
    assert ids == next(csv.reader([stream.readline()], delimiter="\t"))
assert set(cohort.index).issubset(ids) and len(ids) == len(set(ids))
genes = pd.read_csv(file, sep="\t", header=None, skiprows=1, names=["gene", *ids],
                    usecols=["gene", *cohort.index], dtype={"gene": str})
genes = genes.dropna(subset=["gene"]).groupby("gene").median().T
assert {"GZMA", "PRF1"}.issubset(genes.columns)
cyt = genes[["GZMA", "PRF1"]].mean(axis=1, skipna=False)
genes = genes.drop(columns=["GZMA", "PRF1"])
labels = pd.DataFrame({"gene": genes.columns, "feature": [f"P_gene_{i:05d}" for i in range(len(genes.columns))]})
genes.columns = labels.feature
labels.to_csv(inputs / "gene_catalog.csv", index=False)
visits = [f"VOLUME_TUM_BLU_V{k}0" for k in range(1, 5)]
fields = ["CLINICAL-TRIAL-SUBJECT-ID", *visits, "SPHERICITY_T0", "LD_T0", "BPE_5slice_mean_T0"]
mri = pd.read_excel(source / "Multi-feature-MRI-NACT-Data.xlsx", usecols=fields, dtype=str)
mri = mri.set_index(fields[0]).loc[cohort.index].apply(pd.to_numeric, errors="raise")
base = visits[0]
positive = mri.index[mri[base].gt(0) & mri[base].notna()]
clinical = pd.DataFrame(index=positive)
clinical["C_age"] = cohort.Age_at_Screening
clinical["C_hr"] = cohort.HR
clinical["C_her2"] = cohort.HER2
clinical["C_mp"] = cohort.MP
clinical["C_log_volume"] = np.log(mri.loc[positive, base])
clinical["C_longest"] = mri.LD_T0
clinical["C_shape"] = mri.SPHERICITY_T0
clinical["C_background"] = mri.BPE_5slice_mean_T0
clinical["C_cyt"] = cyt
for i, arm in enumerate(sorted(cohort.Arm.dropna().unique())):
    if arm != "Paclitaxel":
        clinical[f"C_arm_{i}"] = cohort.Arm.eq(arm).astype(int)
rows, audit = [], []
for k, col in enumerate(visits[1:], start=1):
    for patient, volume in mri.loc[positive, col].dropna().items():
        if volume < 0:
            raise ValueError("随访体积超出可测范围")
        key = f"ISPY_{patient}_{k}"
        ratio = volume / mri.at[patient, base]
        rows.append({"sample_id": key, "patient": patient, "y": np.log1p(ratio)})
        audit.append({"sample_id": key, "patient": patient, "period": k, "baseline_volume": mri.at[patient, base],
                      "followup_volume": volume})
index = pd.DataFrame(audit).set_index("sample_id")
x = index[["patient"]].join(clinical.join(genes), on="patient").drop(columns="patient")
x["C_period2"] = index.period.eq(2).astype(int)
x["C_period3"] = index.period.eq(3).astype(int)
for period in [2, 3]:
    for variable in ["C_log_volume", "C_hr", "C_her2"]:
        x[f"{variable}_period{period}"] = x[variable] * x[f"C_period{period}"]
x["C_cyt_late"] = x.C_cyt * (index.period - 1) / 2
assert len(genes.columns) > study["maximum_high_features"] and not x.index.duplicated().any()
finish_panel(study, rows, x.reset_index(), audit)
