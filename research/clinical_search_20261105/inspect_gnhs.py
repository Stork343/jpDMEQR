import json
from pathlib import Path

import pandas as pd

root = Path("data/raw/gnhs_20261105/tables")
clinical = pd.read_excel(root / "tableS1" / "tableS1.xlsx", sheet_name="Clinical information")
for column in clinical.columns:
    values = clinical[column]
    print(json.dumps({"column": column, "missing": int(values.isna().sum()), "unique": int(values.nunique()),
                      "examples": values.dropna().astype(str).unique()[:6].tolist()}))
assert clinical.Patient_ID.is_unique
print(json.dumps({"patients": len(clinical), "BMI_complete": clinical[["BMI_BS", "BMI__F2", "BMI__F3"]].notna().sum().to_dict()}))
clinical["canonical_patient"] = "NL" + clinical.Patient_ID.str.removeprefix("NL").astype(int).astype(str)
assert clinical.canonical_patient.is_unique
cohorts = {}
for label in ["tableS2", "tableS3"]:
    matrix = pd.read_excel(root / label / f"{label}.xlsx", sheet_name="B")
    ids = matrix.Sample_ID.str.extract(r"^(NL\d+)(baseline|F2|F3)$")
    assert ids.notna().all().all()
    patient = "NL" + ids[0].str.removeprefix("NL").astype(int).astype(str)
    keys = patient + ids[1]
    assert keys.is_unique
    assert set(patient).issubset(clinical.canonical_patient)
    cohorts[label] = {"patients": set(patient), "samples": set(keys), "proteins": set(matrix.columns[1:])}
    print(json.dumps({"cohort": label, "samples": len(matrix), "patients": patient.nunique(),
                      "visits": ids[1].value_counts().to_dict(), "proteins": len(matrix.columns) - 1}))
print(json.dumps({"overlapping_patients": len(cohorts["tableS2"]["patients"] & cohorts["tableS3"]["patients"]),
                  "overlapping_samples": len(cohorts["tableS2"]["samples"] & cohorts["tableS3"]["samples"]),
                  "common_proteins": len(cohorts["tableS2"]["proteins"] & cohorts["tableS3"]["proteins"])}))
