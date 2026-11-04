import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

config = yaml.safe_load(Path("config/discovery_yeast_stress_20261104.yml").read_text())
root = Path(config["input_directory"])
audit = Path("results/preflight/empirical_selection_20261101/yeast_stress")
root.mkdir(parents=True, exist_ok=True)
audit.mkdir(parents=True, exist_ok=True)
assert not (root / "cohort_counts.csv").exists()
roles = pd.read_csv(Path(config["genomic_directory"]) / "roles.csv")
metadata = pd.read_csv(Path(config["genomic_directory"]) / "metadata.csv")
phenotype = pd.read_csv("data/raw/yeast_1011_20261104/phenoMatrix_35ConditionsNormalizedByYPD.tab.gz", sep="\t", index_col=0)
phenotype.index = phenotype.index.str.removeprefix("SACE_")
assert not phenotype.index.has_duplicates and not roles.strain_id.duplicated().any()
assert set(roles.strain_id) == set(metadata.strain_id)
counts = []
for assay, conditions in config["assays"].items():
    selected = phenotype.loc[roles.strain_id, conditions]
    assert selected.notna().all().all() and np.isfinite(selected.to_numpy()).all()
    assert (selected.to_numpy() >= 0).all()
    rows = selected.rename_axis("strain_id").reset_index().melt(id_vars="strain_id", var_name="condition", value_name="y")
    rows = rows.merge(roles, on="strain_id", validate="many_to_one")
    rows = rows.merge(metadata[["strain_id", "ploidy", "aneuploid"]], on="strain_id", validate="many_to_one")
    rows["sample_id"] = rows.strain_id + "__" + rows.condition
    rows = rows.sort_values(["strain_id", "condition"])
    directory = root / assay
    directory.mkdir()
    rows.loc[rows.role.eq("discovery")].to_csv(directory / "discovery.csv", index=False)
    rows.loc[rows.role.eq("validation")].to_csv(directory / "validation_reserved.csv", index=False)
    table = rows.groupby(["role", "exploration_role"]).agg(strains=("strain_id", "nunique"), records=("sample_id", "size")).reset_index()
    table.insert(0, "assay", assay)
    counts.append(table)
table = pd.concat(counts, ignore_index=True)
table.to_csv(root / "cohort_counts.csv", index=False)
table.to_csv(audit / "cohort_counts.csv", index=False)
(audit / "data_audit.json").write_text(json.dumps({
    "assays": config["assays"], "roles_reused_without_changes": True,
    "formal_validation_used": False, "manuscript_adoption": False,
    "minimum_cluster_size": 2, "maximum_cluster_size": 3,
}, indent=2))
print(table.to_string(index=False))
