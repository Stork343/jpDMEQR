import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

root = Path("data/derived/rice_baap_20261104")
output = Path("results/preflight/empirical_selection_20261101/rice_baap")
output.mkdir(parents=True, exist_ok=True)
config = yaml.safe_load(Path("config/discovery_rice_baap_20261104.yml").read_text(encoding="utf-8"))
assert not (root / "roles.csv").exists()
fam = pd.read_csv(root / "full_genotypes.fam", sep=r"\s+", header=None, dtype=str,
                  names=["FID", "IID", "father", "mother", "sex", "phenotype"])
parts = fam["IID"].str.extract(r"^(\d+)_(\d+)$")
assert parts.notna().all().all() and parts[0].eq(parts[1]).all()
fam["baap_id"] = pd.to_numeric(parts[0], errors="raise")
assert not fam["baap_id"].duplicated().any()
data = pd.read_excel("data/raw/rice_baap_20261104/supplementary_3.xlsx", header=1, na_values=["-"])
data.columns = data.columns.str.strip()
assert data.columns[:2].tolist() == ["BAAP no.", "cultivar name"]
assert not data["BAAP no."].duplicated().any()
conditions = pd.Series(data.columns[2:]).str.extract(r"^Grain cadmium \(ug/kg\)\s+(AWD|CF)\s+year ([12])$")
assert conditions.notna().all().all() and len(conditions) == 4
lookup = pd.DataFrame({"condition": data.columns[2:], "irrigation": conditions[0],
                       "year": 2012 + pd.to_numeric(conditions[1])})
long = data.melt(id_vars=["BAAP no.", "cultivar name"], var_name="condition", value_name="cd_ug_kg")
long = long.merge(lookup, on="condition", validate="many_to_one")
long["cd_ug_kg"] = pd.to_numeric(long["cd_ug_kg"], errors="raise")
observed = long.loc[long["cd_ug_kg"].notna()].copy()
assert observed["cd_ug_kg"].gt(0).all()
matched = observed.merge(fam[["FID", "IID", "baap_id"]], left_on="BAAP no.", right_on="baap_id", validate="many_to_one")
matched = matched.rename(columns={"IID": "accession_id", "cultivar name": "cultivar_name"})
matched["environment"] = matched["irrigation"] + "_" + matched["year"].astype(str)
matched["sample_id"] = matched["accession_id"] + "__" + matched["environment"]
matched["y"] = np.log(matched["cd_ug_kg"])
assert not matched["sample_id"].duplicated().any()
ids = np.sort(matched["accession_id"].unique())
rng = np.random.default_rng(config["role_seed"])
order = rng.permutation(len(ids))
discovery_n = round(len(ids) * config["discovery_fraction"])
training_n = round(discovery_n * config["exploration_training_fraction"])
roles = pd.DataFrame({"accession_id": ids, "role": "validation", "exploration_role": "reserved", "inner_fold": 0})
roles.loc[order[:discovery_n], "role"] = "discovery"
roles.loc[order[:training_n], "exploration_role"] = "training"
roles.loc[order[training_n:discovery_n], "exploration_role"] = "evaluation"
roles.loc[order[:training_n], "inner_fold"] = np.arange(training_n) % config["inner_folds"] + 1
matched = matched.merge(roles, on="accession_id", validate="many_to_one").sort_values(["accession_id", "environment"])
roles.to_csv(root / "roles.csv", index=False)
matched.loc[matched["role"].eq("discovery")].to_csv(root / "discovery.csv", index=False)
matched.loc[matched["role"].eq("validation")].to_csv(root / "validation_reserved.csv", index=False)
training = roles.loc[roles["exploration_role"].eq("training")]
groups = {"training": training["accession_id"].tolist()}
for fold in range(1, config["inner_folds"] + 1):
    groups[f"inner_{fold}"] = training.loc[training["inner_fold"].ne(fold), "accession_id"].tolist()
for name, members in groups.items():
    pd.DataFrame({"FID": "0", "IID": members}).to_csv(root / f"{name}.keep", sep="\t", index=False, header=False)
pd.DataFrame({"FID": "0", "IID": ids}).to_csv(root / "eligible.keep", sep="\t", index=False, header=False)
(root / "training_sets.json").write_text(json.dumps(groups, indent=2), encoding="utf-8")
counts = matched.groupby(["role", "exploration_role"]).agg(records=("sample_id", "size"), accessions=("accession_id", "nunique"))
counts.to_csv(output / "cohort_counts.csv")
cluster_sizes = matched.groupby("accession_id").size()
(output / "data_audit.json").write_text(json.dumps({
    "source_cultivars": len(data), "genotyped_accessions": len(fam), "matched_accessions": len(ids),
    "observed_matched_records": len(matched), "source_missing_records": int(long["cd_ug_kg"].isna().sum()),
    "minimum_cluster_size": int(cluster_sizes.min()), "maximum_cluster_size": int(cluster_sizes.max()),
    "repeated_accessions": int(cluster_sizes.ge(2).sum()),
    "measurement_unit": "published_cultivar_irrigation_year_mean_of_at_least_three_field_replicates",
    "concentration_unit": "microgram_per_kilogram", "response_imputed": False,
    "formal_validation_used": False, "manuscript_adoption": False,
}, indent=2), encoding="utf-8")
print(counts.to_string())
