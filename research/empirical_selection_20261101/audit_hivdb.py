import json
from pathlib import Path

import numpy as np
import pandas as pd

raw = Path("data/raw/hivdb_20261105")
output = Path("results/preflight/empirical_selection_20261101/hivdb_metadata")
output.mkdir(parents=True, exist_ok=True)
reference = pd.read_csv(raw / "PI_DataSet.txt", sep="\t", dtype=str)
full = pd.read_csv(raw / "PI_DataSet.Full.txt", sep="\t", dtype=str)
drugs = ["FPV", "ATV", "IDV", "LPV", "NFV", "SQV", "TPV", "DRV"]
positions = [f"P{i}" for i in range(1, 100)]
assert not reference.SeqID.duplicated().any()
clinical = full.loc[full.Method.eq("PhenoSense") & full.Type.eq("Clinical") & full.SeqID.isin(reference.SeqID)].copy()
known = clinical[positions].fillna(".").apply(lambda column: column.str.fullmatch("[-ACDEFGHIKLMNPQRSTVWY]"))
unique = clinical.loc[~clinical.SeqID.duplicated(keep=False)]
paired = unique[["SeqID"] + drugs].merge(reference[["SeqID"] + drugs], on="SeqID", suffixes=("_full", "_reference"), validate="one_to_one")
records = []
for drug in drugs:
    values = pd.to_numeric(clinical[drug], errors="raise")
    earlier = pd.to_numeric(reference[drug], errors="raise")
    left = pd.to_numeric(paired[f"{drug}_full"], errors="raise")
    right = pd.to_numeric(paired[f"{drug}_reference"], errors="raise")
    observed = left.notna() & right.notna()
    at_100 = observed & left.gt(100) & right.eq(100)
    changed = observed & ~np.isclose(left, right, equal_nan=True)
    records.append({
        "drug": drug, "clinical_numeric_records": int(values.notna().sum()),
        "full_file_maximum": float(values.max()), "reference_file_maximum": float(earlier.max()),
        "matched_full_values_above_100_with_reference_100": int(at_100.sum()),
        "other_matched_numeric_differences": int((changed & ~at_100).sum()),
        "records_equal_to_full_file_maximum": int(values.eq(values.max()).sum()),
    })
pd.DataFrame(records).to_csv(output / "reported_value_audit.csv", index=False)
(output / "status.json").write_text(json.dumps({
    "stage": "metadata_audit", "source_updated": "2026-04-22",
    "reference_isolates": len(reference), "full_records": len(full),
    "clinical_reference_matched_records": len(clinical), "matched_patients": clinical.PtID.nunique(),
    "unambiguous_complete_sequence_patients": clinical.loc[known.all(axis=1), "PtID"].nunique(),
    "unambiguous_complete_B_patients": clinical.loc[known.all(axis=1) & clinical.Subtype.eq("B"), "PtID"].nunique(),
    "instrument_quantification_limit_flags_in_export": False,
    "association_fitted": False, "formal_validation_used": False, "manuscript_adoption": False,
}, indent=2))
print(pd.DataFrame(records).to_string(index=False))
