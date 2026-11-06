from pathlib import Path
import json

import numpy as np
import pandas as pd
import yaml

settings_directory = Path("results/raw/clinical_followup_settings_20261106")
registry = yaml.safe_load(Path("config/clinical_followup_search_20261106.yml").read_text())


def study_config(name):
    return yaml.safe_load((settings_directory / f"{name}.yml").read_text())


def finish_panel(config, rows, features, audit):
    inputs = Path(config["input_directory"])
    report = Path(config["report_directory"])
    roles = pd.read_csv(inputs / "roles.csv", dtype={"patient": str}).set_index("patient")
    rows = pd.DataFrame(rows)
    assert not rows.empty
    rows["patient"] = rows.patient.astype(str)
    assert set(rows.patient).issubset(set(roles.index[roles.role.eq("discovery")]))
    rows["strain_id"] = "PT_" + rows.patient
    rows["role"] = "discovery"
    rows["condition"] = config["response"]
    rows["outer_fold"] = rows.patient.map(roles.outer_fold)
    rows = rows.sort_values(["strain_id", "sample_id"]).reset_index(drop=True)
    x = pd.DataFrame(features).set_index("sample_id").reindex(rows.sample_id)
    assert x.index.is_unique and rows.sample_id.is_unique and np.isfinite(rows.y).all()
    assert rows.y.ge(0).all() and not np.isinf(x.to_numpy(dtype=float)).any()
    patients = rows.groupby("patient").size()
    folds = rows.groupby("outer_fold").patient.nunique().reindex(range(1, 6), fill_value=0)
    target_observed = x[config["inference_target"]].notna().groupby(rows.patient.to_numpy()).mean().mean()
    weight = rows.patient.map(patients).rdiv(1).to_numpy()
    observed = x.notna().mul(weight, axis=0).sum() / len(patients)
    extra = x.columns.str.startswith("P_") & (observed >= config["minimum_observed_fraction"]) & (x.nunique() > 1)
    criteria = {"discovery_patients": len(patients) >= registry["minimum_discovery_patients"],
                "repeated_patients": int(patients.ge(2).sum()) >= registry["minimum_repeated_patients"],
                "fold_patients": int(folds.min()) >= registry["minimum_fold_patients"],
                "distinct_outcomes": rows.y.nunique() >= registry["minimum_distinct_outcomes"],
                "clinical_target_observed": float(target_observed) >= registry["minimum_target_observed_fraction"],
                "additional_features": int(extra.sum()) >= registry["minimum_additional_measurements"]}
    state = {"study_id": config["study_id"], "discovery_patients": len(patients), "records": len(rows),
             "repeated_patients": int(patients.ge(2).sum()), "maximum_records_per_patient": int(patients.max()),
             "distinct_outcomes": rows.y.nunique(), "additional_features_at_registered_frequency": int(extra.sum()),
             "target_observed_fraction": float(target_observed), "criteria": criteria, "passed": all(criteria.values()),
             "reserved_outcomes_used": False, "formal_validation_used": False, "manuscript_adoption": False}
    rows[["sample_id", "strain_id", "y", "condition", "role", "outer_fold"]].to_csv(inputs / "discovery.csv", index=False)
    x.to_csv(inputs / "features.csv.gz")
    pd.DataFrame(audit).to_csv(inputs / "time_audit.csv", index=False)
    folds.rename("patients").to_csv(report / "eligible_by_fold.csv")
    (report / "eligibility.json").write_text(json.dumps(state, indent=2) + "\n")
    print(json.dumps(state), flush=True)
