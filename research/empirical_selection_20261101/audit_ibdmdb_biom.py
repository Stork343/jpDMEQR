import gzip
import json
import shutil
from pathlib import Path

import biom
import numpy as np
import pandas as pd

root = Path("data/raw/ibdmdb_20261101")
compressed = root / "HMP2_metabolomics_w_metadata.biom.gz"
uncompressed = root / "HMP2_metabolomics_w_metadata.biom"
if not uncompressed.exists():
    with gzip.open(compressed, "rb") as source, uncompressed.open("xb") as destination:
        shutil.copyfileobj(source, destination)
table = biom.load_table(str(uncompressed))
assert np.isfinite(table.matrix_data.data).all()
assert (table.matrix_data.data >= 0).all()
output = Path("results/preflight/empirical_selection_20261101/ibdmdb_metadata")
output.mkdir(parents=True, exist_ok=True)
sample_metadata = table.metadata(axis="sample")
feature_metadata = table.metadata(axis="observation")
assert feature_metadata is not None
sample_count = 0 if sample_metadata is None else sum(item is not None for item in sample_metadata)
sample_fields = None if sample_metadata is None else sorted({key for item in sample_metadata if item is not None for key in item})
summary = {
    "features": table.shape[0],
    "samples": table.shape[1],
    "stored_nonzero_values": int(table.matrix_data.nnz),
    "samples_with_embedded_metadata": sample_count,
    "sample_metadata_fields": sample_fields,
    "feature_metadata_fields": sorted({key for item in feature_metadata if item is not None for key in item}),
    "biom_version": biom.__version__,
    "nonnegative_finite_values": True,
    "outcome_associations_fitted": False,
}
(output / "molecular_table.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
annotation = pd.DataFrame(feature_metadata, index=table.ids(axis="observation"))
annotation.to_csv(root / "feature_annotation.csv", index_label="feature_id")
pd.Series(table.ids(), name="sample_id").to_csv(root / "molecular_sample_ids.csv", index=False)
print(json.dumps(summary, indent=2))
for feature_id in table.ids(axis="observation")[:3]:
    print(feature_id, table.metadata(feature_id, axis="observation"))
