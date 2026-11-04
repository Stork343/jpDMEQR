import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from bed_reader import open_bed

raw = Path("data/raw/yeast_1011_20261104")
root = Path("data/derived/yeast_1011_20261104")
audit = Path("results/preflight/empirical_selection_20261101/yeast_1011")
root.mkdir(parents=True, exist_ok=True)
audit.mkdir(parents=True, exist_ok=True)
config = yaml.safe_load(Path("config/discovery_yeast_1011_20261104.yml").read_text())
assert not (root / "roles.csv").exists()


def canonical(values):
    result = pd.Index(values).str.removeprefix("SACE_")
    assert not result.has_duplicates
    return result


phenotype = pd.read_csv(raw / "phenoMatrix_35ConditionsNormalizedByYPD.tab.gz", sep="\t", index_col=0)
phenotype.index = canonical(phenotype.index)
copy = pd.read_csv(raw / "genesMatrix_CopyNumber.tab.gz", sep="\t", index_col=0)
copy.index = canonical(copy.index)
metadata = pd.read_excel(raw / "supplementary_tables.xls", sheet_name="Table S1", header=3)
metadata = metadata.loc[metadata["Standardized name"].notna()].copy()
metadata["strain_id"] = canonical(metadata["Standardized name"])
fam = pd.read_csv(raw / "gwas/1011GWAS_matrix.fam", sep=r"\s+", header=None, dtype=str,
                  names=["FID", "IID", "father", "mother", "sex", "phenotype"])
fam["strain_id"] = canonical(fam["IID"])
assert set(phenotype.index) <= set(fam.strain_id) & set(copy.index) & set(metadata.strain_id)
assert np.isfinite(copy.to_numpy()[~copy.isna().to_numpy()]).all()
assert (copy.fillna(0).to_numpy() >= 0).all()
assert set(metadata["HO deletion"]) == {"yes", "no"}
known_ploidy = metadata["Ploidy"].astype(str).str.fullmatch("[1-5]")
known_aneuploidy = metadata["Aneuploidies"].eq("euploid") | metadata["Aneuploidies"].str.startswith("aneu;", na=False)
eligible = metadata.loc[metadata["HO deletion"].eq("no") & known_ploidy & known_aneuploidy].copy()
eligible = eligible.loc[eligible.strain_id.isin(phenotype.index)].sort_values("strain_id")
eligible["ploidy"] = pd.to_numeric(eligible["Ploidy"], errors="raise")
eligible["aneuploid"] = eligible["Aneuploidies"].ne("euploid").astype(int)
ids = eligible.strain_id.to_numpy()
values = phenotype.loc[ids, config["conditions"]]
assert values.notna().all().all() and np.isfinite(values.to_numpy()).all()
assert (values.to_numpy() >= 0).all()
roles = pd.DataFrame({"strain_id": ids, "role": "validation", "exploration_role": "reserved", "inner_fold": 0})
rng = np.random.default_rng(config["role_seed"])
order = rng.permutation(len(ids))
discovery_n = round(len(ids) * config["discovery_fraction"])
training_n = round(discovery_n * config["exploration_training_fraction"])
roles.loc[order[:discovery_n], "role"] = "discovery"
roles.loc[order[:training_n], "exploration_role"] = "training"
roles.loc[order[training_n:discovery_n], "exploration_role"] = "evaluation"
roles.loc[order[:training_n], "inner_fold"] = np.arange(training_n) % config["inner_folds"] + 1
long = values.rename_axis("strain_id").reset_index().melt(id_vars="strain_id", var_name="condition", value_name="y")
long["high_salt"] = long.condition.eq("YPDNACL15M").astype(int)
long["sample_id"] = long.strain_id + "__" + long.condition
long = long.merge(eligible[["strain_id", "ploidy", "aneuploid"]], on="strain_id", validate="many_to_one")
long = long.merge(roles, on="strain_id", validate="many_to_one").sort_values(["strain_id", "condition"])
assert long.groupby("strain_id").size().eq(2).all()
roles.to_csv(root / "roles.csv", index=False)
long.loc[long.role.eq("discovery")].to_csv(root / "discovery.csv", index=False)
long.loc[long.role.eq("validation")].to_csv(root / "validation_reserved.csv", index=False)
copy.loc[ids].rename_axis("strain_id").to_csv(root / "copy_number.csv.gz")
fam[["FID", "IID", "strain_id"]].to_csv(root / "genotype_identifiers.csv", index=False)
eligible[["strain_id", "Isolate name", "Standardized name", "Clades", "ploidy", "aneuploid"]].to_csv(root / "metadata.csv", index=False)
training = roles.loc[roles.exploration_role.eq("training")]
groups = {"training": training.strain_id.tolist()}
for fold in range(1, config["inner_folds"] + 1):
    groups[f"inner_{fold}"] = training.loc[training.inner_fold.ne(fold), "strain_id"].tolist()
for name, members in {**groups, "eligible": ids.tolist()}.items():
    fam.loc[fam.strain_id.isin(members), ["FID", "IID"]].to_csv(root / f"{name}.keep", sep="\t", index=False, header=False)
(root / "training_sets.json").write_text(json.dumps(groups, indent=2))

bed = open_bed(raw / "gwas/1011GWAS_matrix.bed", count_A1=True)
bim = pd.read_csv(raw / "gwas/1011GWAS_matrix.bim", sep=r"\s+", header=None)
sample_index = np.linspace(0, len(fam) - 1, 12, dtype=int)
variant_index = np.linspace(0, bim[0].between(1, 16).sum() - 1, 64, dtype=int)
reference = pd.DataFrame(bed.read(index=(sample_index, variant_index)), columns=bim.iloc[variant_index, 1])
reference.insert(0, "IID", fam.iloc[sample_index].IID.to_numpy())
reference.to_csv(root / "genotype_reader_reference.csv", index=False)
counts = long.groupby(["role", "exploration_role"]).agg(strains=("strain_id", "nunique"), records=("sample_id", "size"))
counts.to_csv(audit / "cohort_counts.csv")
eligible.groupby(["ploidy", "aneuploid"]).size().rename("strains").to_csv(audit / "cohort_metadata.csv")
prior = pd.read_excel(raw / "supplementary_tables.xls", sheet_name="Table S15", header=2)
prior.loc[prior.Condition.isin(config["conditions"])].to_csv(audit / "original_salt_findings.csv", index=False)
(audit / "data_audit.json").write_text(json.dumps({
    "source_phenotyped_strains": len(phenotype), "source_genotyped_strains": len(fam),
    "eligible_strains": len(ids), "records": len(long), "cluster_size": 2,
    "measurement": "reported_strain_condition_growth_ratio", "raw_colony_replicates_available": False,
    "copy_number_columns": len(copy.columns), "copy_number_scale": "estimated_per_haploid_genome",
    "formal_validation_used": False, "manuscript_adoption": False,
    "dependencies": {"xlrd": "2.0.2", "bed_reader": "1.1.0"},
}, indent=2))
print(counts.to_string())
