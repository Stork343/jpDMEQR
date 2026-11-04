import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import yaml

raw = Path("data/raw/arabidopsis_20261105")
root = Path("data/derived/arabidopsis_20261105")
audit = Path("results/preflight/empirical_selection_20261101/arabidopsis")
root.mkdir(parents=True, exist_ok=True)
audit.mkdir(parents=True, exist_ok=True)
assert not (root / "roles.csv").exists()
config = yaml.safe_load(Path("config/discovery_arabidopsis_20261105.yml").read_text())
with h5py.File(raw / "1001_SNP_MATRIX/imputed_snps_binary.hdf5", "r") as source:
    genotype_ids = source["accessions"].asstr()[:]
    assert len(set(genotype_ids)) == len(genotype_ids)
tables = [pd.read_csv(raw / f"{condition}.csv", dtype={"accession_id": str}) for condition in config["conditions"]]
paired = sorted(set(tables[0].accession_id) & set(tables[1].accession_id) & set(genotype_ids))
rows = pd.concat(tables, ignore_index=True)
assert not rows[["phenotype_name", "obs_unit_id"]].duplicated().any()
rows = rows.loc[rows.accession_id.isin(paired)].copy()
rows = rows.rename(columns={"phenotype_name": "condition", "phenotype_value": "y"})
rows["high_temperature"] = rows.condition.eq("FT16").astype(int)
assert rows.y.notna().all() and np.isfinite(rows.y).all() and rows.y.gt(0).all()
counts_by_condition = rows.groupby(["accession_id", "condition"]).size().unstack()
assert counts_by_condition.FT10.eq(counts_by_condition.FT16).all()
assert counts_by_condition.to_numpy().min() == 1 and counts_by_condition.to_numpy().max() <= 2
repeated = rows.loc[rows.accession_id.isin(counts_by_condition.index[counts_by_condition.FT10.gt(1)])]
repeated.to_csv(root / "repeated_source_units.csv", index=False)
rng = np.random.default_rng(config["role_seed"])
order = rng.permutation(len(paired))
discovery_n = round(len(paired) * config["discovery_fraction"])
training_n = round(discovery_n * config["exploration_training_fraction"])
roles = pd.DataFrame({"accession_id": paired, "role": "validation", "exploration_role": "reserved", "inner_fold": 0})
roles.loc[order[:discovery_n], "role"] = "discovery"
roles.loc[order[:training_n], "exploration_role"] = "training"
roles.loc[order[training_n:discovery_n], "exploration_role"] = "evaluation"
roles.loc[order[:training_n], "inner_fold"] = np.arange(training_n) % config["inner_folds"] + 1
rows = rows.merge(roles, on="accession_id", validate="many_to_one").sort_values(["accession_id", "condition", "obs_unit_id"])
rows["sample_id"] = rows.accession_id + "__" + rows.condition + "__" + rows.obs_unit_id.astype(str)
roles.to_csv(root / "roles.csv", index=False)
rows.loc[rows.role.eq("discovery")].to_csv(root / "discovery.csv", index=False)
rows.loc[rows.role.eq("validation")].to_csv(root / "validation_reserved.csv", index=False)
training = roles.loc[roles.exploration_role.eq("training")]
groups = {"training": training.accession_id.tolist()}
for fold in range(1, config["inner_folds"] + 1):
    groups[f"inner_{fold}"] = training.loc[training.inner_fold.ne(fold), "accession_id"].tolist()
for name, members in {**groups, "eligible": paired}.items():
    pd.DataFrame({"FID": "0", "IID": members}).to_csv(root / f"{name}.keep", sep="\t", index=False, header=False)
(root / "training_sets.json").write_text(json.dumps(groups, indent=2))
summary = rows.groupby(["role", "exploration_role"]).agg(accessions=("accession_id", "nunique"), records=("sample_id", "size"))
summary.to_csv(audit / "cohort_counts.csv")
(audit / "data_audit.json").write_text(json.dumps({
    "source_FT10_records": len(tables[0]), "source_FT16_records": len(tables[1]),
    "source_genotyped_accessions": len(genotype_ids), "matched_accessions": len(paired),
    "records": len(rows), "measurement": "published_mean_of_flowering_time_replicates",
    "original_protocol_replicates": 4, "genotype_source": "author_imputed_binary_SNP_matrix",
    "source_units_retained_for_accession_8424": True, "maximum_cluster_size": int(rows.groupby("accession_id").size().max()),
    "response_imputed": False, "formal_validation_used": False, "manuscript_adoption": False,
}, indent=2))
print(summary.to_string())
