import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd
import yaml

config = yaml.safe_load(Path("config/discovery_yeast_proteome_20261105.yml").read_text())
root = Path(config["input_directory"])
raw = Path("data/raw/yeast_1011_20261104")
audit = Path("results/preflight/empirical_selection_20261101/yeast_proteome")
assert not (root / "cohort_counts.csv").exists()
roles = pd.read_csv(root / "roles.csv")
phenotype = pd.read_csv(raw / "phenoMatrix_35ConditionsNormalizedByYPD.tab.gz", sep="\t", index_col=0)
phenotype.index = phenotype.index.str.removeprefix("SACE_")
metadata = pd.read_csv("data/derived/yeast_1011_20261104/metadata.csv")
with ZipFile(raw / "Teyssonniere_2024_supplements.zip") as source:
    assert source.testzip() is None
    cultivation = pd.read_excel(source.open("pnas.2319211121.sd01.xlsx"))
    published = pd.read_excel(source.open("pnas.2319211121.sd02.xlsx"), index_col=0)
    source.extract("pnas.2319211121.sd01.xlsx", raw / "Teyssonniere_2024")
    source.extract("pnas.2319211121.sd02.xlsx", raw / "Teyssonniere_2024")
cultivation["strain_id"] = cultivation.Standardized_name.str.removeprefix("SACE_")
published.columns = published.columns.str.replace(".", "_", regex=False).str.removeprefix("SACE_")
protein = pd.read_csv(root / "protein_abundance.csv.gz", index_col=0)
expected = published.loc[protein.index, protein.columns]
assert np.array_equal(expected.isna(), protein.isna())
relative_error = np.nanmax(np.abs(expected.to_numpy() - protein.to_numpy()) / expected.to_numpy())
assert relative_error < 1e-12
information = metadata[["strain_id", "ploidy", "aneuploid"]].merge(
    cultivation[["strain_id", "OD_proteomic"]], on="strain_id", validate="one_to_one")
information = information.rename(columns={"OD_proteomic": "harvest_od"})
assert set(roles.strain_id) <= set(information.strain_id)
matched = information.loc[information.strain_id.isin(roles.strain_id)]
assert matched.harvest_od.notna().all() and matched.harvest_od.gt(0).all()
counts = []
for assay_name, conditions in config["assays"].items():
    values = phenotype.loc[roles.strain_id, conditions]
    assert values.notna().all().all() and np.isfinite(values.to_numpy()).all()
    rows = values.rename_axis("strain_id").reset_index().melt(id_vars="strain_id", var_name="condition", value_name="y")
    rows = rows.merge(roles, on="strain_id", validate="many_to_one")
    rows = rows.merge(information, on="strain_id", validate="many_to_one")
    rows["sample_id"] = rows.strain_id + "__" + rows.condition
    rows = rows.sort_values(["strain_id", "condition"])
    directory = root / assay_name
    directory.mkdir()
    rows.loc[rows.role.eq("discovery")].to_csv(directory / "discovery.csv", index=False)
    rows.loc[rows.role.eq("validation")].to_csv(directory / "validation_reserved.csv", index=False)
    table = rows.groupby(["role", "exploration_role"]).agg(strains=("strain_id", "nunique"), records=("sample_id", "size")).reset_index()
    table.insert(0, "assay", assay_name)
    counts.append(table)
table = pd.concat(counts, ignore_index=True)
table.to_csv(root / "cohort_counts.csv", index=False)
table.to_csv(audit / "analysis_cohort_counts.csv", index=False)
(audit / "source_checks.json").write_text(json.dumps({
    "status": "passed", "published_matrix_relative_error": float(relative_error),
    "cultivation_OD_observed": len(matched), "cultivation_OD_missing": 0,
    "formal_validation_used": False, "manuscript_adoption": False,
}, indent=2))
print(table.to_string(index=False))
