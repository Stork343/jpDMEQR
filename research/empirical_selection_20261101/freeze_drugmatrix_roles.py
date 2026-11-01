import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

config = yaml.safe_load(Path("config/discovery_drugmatrix_20261101.yml").read_text(encoding="utf-8"))
data = pd.read_csv("data/raw/drugmatrix_20261101/matched_liver_animals.csv")
data = data.loc[
    data["STUDY_LAB"].eq("SRI International")
    & data["SEX_array"].eq("M")
    & data["ROUTE_clinical"].eq("Oral Gavage")
    & data["TREATMENT"].eq("Chemical")
].copy()
assert data["DOSE_array"].str.endswith(" mg/kg").all()
array_dose = pd.to_numeric(data["DOSE_array"].str.removesuffix(" mg/kg"), errors="raise")
assert np.allclose(array_dose, data["DOSE_clinical"], rtol=1e-8, atol=1e-8)
assert not data["ANIMAL_ID"].duplicated().any()
compounds = np.sort(data["TEST_ARTICLE"].unique())
rng = np.random.default_rng(config["role_seed"])
order = rng.permutation(len(compounds))
discovery_n = round(len(compounds) * config["discovery_fraction"])
training_n = round(discovery_n * config["exploration_training_fraction"])
roles = pd.DataFrame({"compound_id": compounds, "role": "validation", "exploration_role": "reserved", "inner_fold": 0})
roles.loc[order[:discovery_n], "role"] = "discovery"
roles.loc[order[:training_n], "exploration_role"] = "training"
roles.loc[order[training_n:discovery_n], "exploration_role"] = "evaluation"
roles.loc[order[:training_n], "inner_fold"] = np.arange(training_n) % config["inner_folds"] + 1
directory = Path("data/derived/drugmatrix_20261101")
directory.mkdir(parents=True, exist_ok=True)
assert not (directory / "roles.csv").exists()
roles.to_csv(directory / "roles.csv", index=False)
data = data.merge(roles, left_on="TEST_ARTICLE", right_on="compound_id", validate="many_to_one")
data.to_csv(directory / "cohort_metadata.csv", index=False)
output = Path("results/preflight/empirical_selection_20261101/drugmatrix_metadata")
counts = data.groupby(["role", "exploration_role"]).agg(animals=("ANIMAL_ID", "size"), compounds=("compound_id", "nunique"))
counts.to_csv(output / "frozen_cohort.csv")
(output / "role_freeze.json").write_text(json.dumps({
    "role_seed": config["role_seed"], "compounds": len(compounds), "animals": len(data),
    "dose_metadata_consistent": True, "outcome_associations_fitted": False,
    "formal_validation_used": False, "manuscript_adoption": False,
}, indent=2), encoding="utf-8")
print(counts.to_string())
