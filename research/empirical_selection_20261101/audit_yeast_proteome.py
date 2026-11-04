import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd

raw = Path("data/raw/yeast_1011_20261104")
root = Path("data/derived/yeast_proteome_20261104")
audit = Path("results/preflight/empirical_selection_20261101/yeast_proteome")
root.mkdir(parents=True, exist_ok=True)
audit.mkdir(parents=True, exist_ok=True)
archive = raw / "Teyssonniere_data_upload.zip"
member = "data_upload/DIA-NN_1.8/SCmedia_MBR_CommonPeptide_Approach/211025_ProteomicsData_SCmediaMBR_DIA-NN-1.8_genes_ORF.tsv"
with ZipFile(archive) as source:
    assert source.testzip() is None
    protein = pd.read_csv(source.open(member), sep="\t", index_col="Protein.Group")
    annotation = pd.read_csv(source.open("data_upload/CompleteGeneAnnot_strand.csv"))
    source.extract(member, raw / "Teyssonniere_2024")
    source.extract("data_upload/CompleteGeneAnnot_strand.csv", raw / "Teyssonniere_2024")
protein.columns = protein.columns.str.replace(".", "_", regex=False).str.removeprefix("SACE_")
assert protein.shape == (630, 942)
assert not protein.columns.has_duplicates and not protein.index.has_duplicates
observed = protein.to_numpy()[protein.notna().to_numpy()]
assert np.isfinite(observed).all() and (observed > 0).all()
roles = pd.read_csv("data/derived/yeast_1011_20261104/roles.csv")
matched = roles.loc[roles.strain_id.isin(protein.columns)].copy()
matrix = protein.loc[:, matched.strain_id]
matrix.to_csv(root / "protein_abundance.csv.gz", index_label="protein_id")
matched.to_csv(root / "roles.csv", index=False)
counts = matched.groupby(["role", "exploration_role"]).size().rename("strains").reset_index()
counts.to_csv(audit / "cohort_counts.csv", index=False)
missing = matrix.isna().mean(axis=0)
(audit / "data_audit.json").write_text(json.dumps({
    "source_doi": "10.1073/pnas.2319211121", "data_doi": "10.5281/zenodo.10567083",
    "archive_member": member, "source_proteins": len(protein), "source_strains": protein.shape[1],
    "matched_strains": len(matched), "missing_fraction": float(matrix.isna().to_numpy().mean()),
    "maximum_strain_missing_fraction": float(missing.max()),
    "source_processing": "DIA_NN_maxLFQ_with_published_precursor_median_batch_adjustment",
    "source_imputed": False, "source_cross_layer_quantile_normalized": False,
    "roles_inherited": True, "formal_validation_used": False, "association_fitted": False,
    "annotation_columns": annotation.columns.tolist(), "manuscript_adoption": False,
}, indent=2))
print(counts.to_string(index=False))
print("Missing fraction", matrix.isna().to_numpy().mean())
print("Annotation columns", annotation.columns.tolist())
