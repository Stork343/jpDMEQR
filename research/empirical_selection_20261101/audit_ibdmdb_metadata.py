import json
from pathlib import Path

import pandas as pd

root = Path("data/raw/ibdmdb_20261101")
output = Path("results/preflight/empirical_selection_20261101/ibdmdb_metadata")
output.mkdir(parents=True, exist_ok=True)
data = pd.read_csv(root / "hmp2_metadata_2018-08-20.csv", low_memory=False)
required = ["External ID", "Participant ID", "data_type", "week_num", "diagnosis", "fecalcal"]
assert set(required).issubset(data.columns)
data["fecalcal"] = pd.to_numeric(data["fecalcal"], errors="raise")
counts = data.groupby("data_type", dropna=False).agg(
    rows=("External ID", "size"),
    samples=("External ID", "nunique"),
    participants=("Participant ID", "nunique"),
    calprotectin_records=("fecalcal", "count"),
)
counts.to_csv(output / "assay_counts.csv")
selected = data.loc[data["data_type"].str.contains("metabol", case=False, na=False)].copy()
available = selected.dropna(subset=["fecalcal"])
repeats = available.groupby("Participant ID")["External ID"].nunique()
values = available["fecalcal"]
summary = {
    "stage": "metadata_audit",
    "metabolomics_records": len(selected),
    "matched_calprotectin_records": len(available),
    "matched_participants": int(available["Participant ID"].nunique()),
    "participants_with_repeats": int((repeats >= 2).sum()),
    "duplicate_metabolomics_sample_ids": int(selected["External ID"].duplicated().sum()),
    "nonpositive_calprotectin": int((values <= 0).sum()),
    "distinct_calprotectin_values": int(values.nunique()),
    "calprotectin_minimum": float(values.min()) if len(values) else None,
    "calprotectin_maximum": float(values.max()) if len(values) else None,
    "largest_tied_fraction": float(values.value_counts().max() / len(values)) if len(values) else None,
    "maximum_measurements": int(repeats.max()) if len(repeats) else None,
    "outcome_associations_fitted": False,
    "manuscript_adoption": False,
}
(output / "audit.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(counts.to_string())
print(json.dumps(summary, indent=2))
