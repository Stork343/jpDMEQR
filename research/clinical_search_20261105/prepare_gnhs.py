import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

parser = argparse.ArgumentParser()
parser.add_argument("config", type=Path)
args = parser.parse_args()
config = yaml.safe_load(args.config.read_text())
source = Path(config["source_directory"]) / "tables"
output = Path(config["input_directory"])
assert not output.exists()
clinical = pd.read_excel(source / "tableS1" / "tableS1.xlsx", sheet_name="Clinical information")
clinical["strain_id"] = "NL" + clinical.Patient_ID.str.removeprefix("NL").astype(int).astype(str)
assert clinical.strain_id.is_unique
clinical = clinical.set_index("strain_id")
matrices = {}
people = {}
for label, role in [("tableS2", "discovery"), ("tableS3", "validation")]:
    frame = pd.read_excel(source / label / f"{label}.xlsx", sheet_name="B")
    parts = frame.Sample_ID.str.extract(r"^(NL\d+)(baseline|F2|F3)$")
    assert parts.notna().all().all()
    patient = "NL" + parts[0].str.removeprefix("NL").astype(int).astype(str)
    frame["strain_id"] = patient
    frame["visit"] = parts[1]
    assert not frame[["strain_id", "visit"]].duplicated().any()
    matrices[role] = frame
    people[role] = set(patient)
assert not people["discovery"] & people["validation"]
proteins = sorted((set(matrices["discovery"].columns) & set(matrices["validation"].columns)) -
                  {"Sample_ID", "strain_id", "visit"})
assert len(proteins) == 413
rows = []
features = []
exclusions = []
bmi = {"baseline": "BMI_BS", "F2": "BMI__F2", "F3": "BMI__F3"}
age = {"baseline": "Age_BS", "F2": "Age__F2", "F3": "Age__F3"}
intervals = [("baseline", "F2"), ("baseline", "F3")] if config.get("prediction_origin") == "baseline" else [("baseline", "F2"), ("F2", "F3")]
for role, frame in matrices.items():
    for start, end in intervals:
        current = frame.loc[frame.visit.eq(start)].copy()
        phenotypes = clinical.loc[current.strain_id].copy()
        eligible = phenotypes[bmi[start]].notna() & phenotypes[bmi[end]].notna()
        exclusions.append({"role": role, "origin": start, "available_protein_records": len(current),
                           "missing_prior_or_future_bmi": int((~eligible).sum())})
        current = current.loc[eligible.to_numpy()]
        phenotypes = phenotypes.loc[eligible]
        prior = phenotypes[bmi[start]].to_numpy(dtype=float)
        outcome = phenotypes[bmi[end]].to_numpy(dtype=float)
        start_age = phenotypes[age[start]].to_numpy(dtype=float)
        horizon = phenotypes[age[end]].to_numpy(dtype=float) - start_age
        assert np.isfinite(prior).all() and np.isfinite(outcome).all() and np.isfinite(horizon).all()
        assert ((prior > 5) & (prior < 100) & (outcome > 5) & (outcome < 100)).all()
        assert (horizon > 0).all()
        sample_id = current.strain_id + "_" + start + "_to_" + end
        rows.append(pd.DataFrame({"sample_id": sample_id, "strain_id": current.strain_id,
                                 "y": outcome, "condition": "future_bmi", "role": role}))
        C = pd.DataFrame({"sample_id": sample_id, "C_prior_bmi": prior, "C_prior_age": start_age,
                          "C_sex_code": phenotypes.Sex.to_numpy(), "C_horizon": horizon,
                          "C_origin_F2": float(start == "F2"), "C_bmi_horizon": prior * horizon,
                          "C_age_horizon": start_age * horizon}).set_index("sample_id")
        if config.get("prediction_origin") == "baseline":
            C["C_target_F3"] = float(end == "F3")
        suffix = "BS" if start == "baseline" else "F2"
        for name in ["Waist", "Glucose", "TG", "HDL", "SBP/DBP"]:
            values = phenotypes[f"{name}_{suffix}"].replace(999, np.nan)
            assert values.dropna().isin([0, 1]).all()
            C["C_prior_" + name.replace("/", "_")] = values.to_numpy()
        P = current[proteins].copy()
        assert np.isfinite(P.to_numpy()[P.notna().to_numpy()]).all()
        P.columns = ["P_" + column for column in proteins]
        P.index = sample_id
        features.append(C.join(P, validate="one_to_one"))
rows = pd.concat(rows, ignore_index=True).sort_values(["strain_id", "sample_id"])
features = pd.concat(features).loc[rows.sample_id]
assert rows.sample_id.is_unique and features.index.is_unique
if config.get("prediction_origin") == "baseline":
    fixed = features[[column for column in features.columns if column.startswith("P_") or column.startswith("C_prior_") or column == "C_sex_code"]].copy()
    fixed["strain_id"] = rows.strain_id.to_numpy()
    assert fixed.groupby("strain_id").nunique(dropna=False).le(1).all().all()
participants = rows[["strain_id", "role"]].drop_duplicates().sort_values("strain_id").reset_index(drop=True)
assert participants.strain_id.is_unique
participants["outer_fold"] = 0
participants["sex_code"] = clinical.loc[participants.strain_id, "Sex"].to_numpy()
rng = np.random.default_rng(config["role_seed"])
for sex, indices in participants.loc[participants.role.eq("discovery")].groupby("sex_code").groups.items():
    selected = rng.permutation(np.asarray(indices))
    participants.loc[selected, "outer_fold"] = np.arange(len(selected)) % config["outer_folds"] + 1
rows = rows.merge(participants[["strain_id", "outer_fold"]], on="strain_id", validate="many_to_one")
output.mkdir(parents=True)
participants.to_csv(output / "roles.csv", index=False)
features.to_csv(output / "features.csv.gz", index_label="sample_id")
pd.DataFrame({"protein": proteins}).to_csv(output / "protein_catalog.csv", index=False)
for role in ["discovery", "validation"]:
    rows.loc[rows.role.eq(role)].to_csv(output / f"{role}.csv", index=False)
report = Path(config["report_directory"])
report.mkdir(parents=True, exist_ok=True)
counts = rows.groupby("role").agg(patients=("strain_id", "nunique"), records=("sample_id", "size"))
counts.to_csv(report / "cohort_counts.csv")
pd.DataFrame(exclusions).to_csv(report / "cohort_exclusions.csv", index=False)
sizes = rows.groupby("strain_id").size()
status = {"status": "prepared", "source_doi": config["source_doi"], "patients": rows.strain_id.nunique(),
          "records": len(rows), "patients_with_repeats": int(sizes.ge(2).sum()), "maximum_cluster_size": int(sizes.max()),
          "protein_candidates": len(proteins), "role_assignment": config["role_assignment"],
          "prediction_origin": config.get("prediction_origin", "adjacent_visits"),
          "patient_overlap": 0, "outcome_imputed": False, "future_assay_used_as_predictor": False,
          "formal_validation_used": False, "manuscript_adoption": False}
(report / "cohort_status.json").write_text(json.dumps(status, indent=2) + "\n")
manifest = [json.loads((Path(config["source_directory"]) / f"tableS{i}.source.json").read_text()) for i in [1, 2, 3]]
(report / "data_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
print(json.dumps(status, indent=2))
print(counts.to_string())
