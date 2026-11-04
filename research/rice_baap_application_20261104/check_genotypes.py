import json
from pathlib import Path

import numpy as np
import pandas as pd
from bed_reader import open_bed

root = Path("data/derived/rice_baap_20261104")
output = Path("results/scratch/rice_baap_20261104")
bed = open_bed(root / "training_panel.bed", count_A1=True)
sample_index = np.linspace(0, bed.iid_count - 1, 12, dtype=int)
variant_index = np.linspace(0, bed.sid_count - 1, 64, dtype=int)
values = bed.read(index=(sample_index, variant_index), dtype="float64", num_threads=1)
assert np.isin(values, [0, 1, 2]).all()
pd.DataFrame(values, index=bed.iid[sample_index], columns=bed.sid[variant_index]).to_csv(
    output / "bed_reader_reference.csv", index_label="accession_id"
)
print(json.dumps({"samples": bed.iid_count, "markers": bed.sid_count,
                  "reference_genotypes": values.size, "counts": "A1"}))
