import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

config = yaml.safe_load(Path("config/application_gnhs_future_bmi_20261105.yml").read_text())
source = Path(config["source_directory"]) / "tables"
directory = Path(config["input_directory"])
clinical = pd.read_excel(source / "tableS1" / "tableS1.xlsx", sheet_name="Clinical information")
clinical["strain_id"] = "NL" + clinical.Patient_ID.str.removeprefix("NL").astype(int).astype(str)
bmi = clinical.melt(id_vars="strain_id", value_vars=["BMI_BS", "BMI__F2", "BMI__F3"], value_name="bmi")
bmi["visit"] = bmi.variable.map({"BMI_BS": "baseline", "BMI__F2": "F2", "BMI__F3": "F3"})
bmi = bmi.set_index(["strain_id", "visit"]).bmi
age = clinical.melt(id_vars="strain_id", value_vars=["Age_BS", "Age__F2", "Age__F3"], value_name="age")
age["visit"] = age.variable.map({"Age_BS": "baseline", "Age__F2": "F2", "Age__F3": "F3"})
age = age.set_index(["strain_id", "visit"]).age
features = pd.read_csv(directory / "features.csv.gz").set_index("sample_id")
catalog = pd.read_csv(directory / "protein_catalog.csv")
records = []
people = {}
for label, role in [("tableS2", "discovery"), ("tableS3", "validation")]:
    rows = pd.read_csv(directory / f"{role}.csv")
    parts = rows.sample_id.str.extract(r"^(NL\d+)_(baseline|F2)_to_(F2|F3)$")
    assert parts.notna().all().all() and parts[0].equals(rows.strain_id)
    pairs = list(zip(parts[1], parts[2]))
    assert set(pairs).issubset({("baseline", "F2"), ("F2", "F3")})
    start_key = pd.MultiIndex.from_arrays([rows.strain_id, parts[1]])
    end_key = pd.MultiIndex.from_arrays([rows.strain_id, parts[2]])
    expected_y = bmi.loc[end_key].to_numpy()
    expected_prior = bmi.loc[start_key].to_numpy()
    expected_age = age.loc[start_key].to_numpy()
    expected_horizon = age.loc[end_key].to_numpy() - expected_age
    stored = features.loc[rows.sample_id]
    clinical_error = np.max(np.abs(np.column_stack([
        expected_y - rows.y, expected_prior - stored.C_prior_bmi.to_numpy(),
        expected_age - stored.C_prior_age.to_numpy(), expected_horizon - stored.C_horizon.to_numpy()])))
    assert np.isfinite(clinical_error) and clinical_error < 1e-10
    matrix = pd.read_excel(source / label / f"{label}.xlsx", sheet_name="B")
    ids = matrix.Sample_ID.str.extract(r"^(NL\d+)(baseline|F2|F3)$")
    assert ids.notna().all().all()
    matrix["patient"] = "NL" + ids[0].str.removeprefix("NL").astype(int).astype(str)
    matrix["visit"] = ids[1]
    matrix = matrix.set_index(["patient", "visit"])
    expected_protein = matrix.loc[start_key, catalog.protein].to_numpy(dtype=float)
    stored_protein = stored[["P_" + name for name in catalog.protein]].to_numpy(dtype=float)
    assert np.array_equal(np.isnan(expected_protein), np.isnan(stored_protein))
    protein_error = float(np.nanmax(np.abs(expected_protein - stored_protein)))
    assert protein_error < 1e-10
    assert rows.groupby("strain_id").size().max() <= 2
    people[role] = set(rows.strain_id)
    records.append({"role": role, "patients": rows.strain_id.nunique(), "records": len(rows),
                    "clinical_alignment_error": float(clinical_error), "protein_alignment_error": protein_error,
                    "predictor_visit_precedes_response": bool((expected_horizon > 0).all())})
assert not people["discovery"] & people["validation"]
for outer in range(1, config["outer_folds"] + 1):
    fold = Path(config["output_root"]) / f"outer_{outer}"
    train = pd.read_csv(fold / "training.csv")
    test = pd.read_csv(fold / "evaluation.csv")
    assert not set(train.strain_id) & set(test.strain_id)
    assert set(train.strain_id) | set(test.strain_id) == people["discovery"]
    assert not (set(train.strain_id) | set(test.strain_id)) & people["validation"]
output = Path(config["report_directory"])
pd.DataFrame(records).to_csv(output / "temporal_alignment_checks.csv", index=False)
status = {"status": "passed", "patient_overlap": 0, "folds_checked": config["outer_folds"],
          "records_checked": sum(record["records"] for record in records), "formal_validation_used": False}
(output / "temporal_alignment_status.json").write_text(json.dumps(status, indent=2) + "\n")
print(json.dumps(records, indent=2))
