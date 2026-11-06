import csv
import gzip
import json
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd

from panel import study_config

study = study_config("ispy2_volume")
inputs = Path(study["input_directory"])
source = Path(study["source_directory"])
roles = pd.read_csv(inputs / "roles.csv", dtype={"patient": str}).set_index("patient")
discovery = set(roles.index[roles.role.eq("discovery")])
with gzip.open(source / "expression.tsv.gz", "rt") as stream:
    records = csv.reader(stream, delimiter="\t")
    patients = next(records)
    indices = [(i+1, k) for i,k in enumerate(patients) if k in discovery]
    programs = {"PRF1": [], "GZMA": []}
    for record in records:
        if record[0] in programs:
            programs[record[0]].append([float(record[i]) for i, _ in indices])
assert all(programs.values())
activity = (np.median(programs["PRF1"], axis=0) + np.median(programs["GZMA"], axis=0)) / 2
activity = dict(zip([k for _, k in indices], activity))
book = openpyxl.load_workbook(source / "Multi-feature-MRI-NACT-Data.xlsx", read_only=True, data_only=True)
sheet, = book.worksheets
records = sheet.values
fields = list(next(records))
expected = []
for row in records:
    patient = str(row[fields.index("CLINICAL-TRIAL-SUBJECT-ID")])
    if patient not in discovery:
        continue
    baseline = row[fields.index("VOLUME_TUM_BLU_V10")]
    if baseline is None or baseline <= 0:
        continue
    for period in (1, 2, 3):
        volume = row[fields.index(f"VOLUME_TUM_BLU_V{period+1}0")]
        if volume is not None:
            ratio = volume/baseline
            expected.append({"sample_id": f"ISPY_{patient}_{period}", "reference_y": np.log1p(ratio),
                             "ratio": ratio, "published_change": row[fields.index(f"FTV_pch_T0_T{period}")],
                             "score": activity[patient]*(period-1)/2, "patient": patient})
book.close()
reference = pd.DataFrame(expected)
rows = pd.read_csv(inputs / "discovery.csv")
x = pd.read_csv(inputs / "features.csv.gz", usecols=["sample_id", "C_cyt_late"])
paired = rows.merge(reference, on="sample_id", validate="one_to_one").merge(x, on="sample_id", validate="one_to_one")
assert len(paired) == len(rows) == len(reference)
np.testing.assert_allclose(paired.y, paired.reference_y, atol=1e-12)
np.testing.assert_allclose(paired.C_cyt_late, paired.score, atol=1e-12)
changes = paired.dropna(subset=["published_change"])
errors = {unit: float((changes.published_change-factor*(changes.ratio-1)).abs().max())
          for unit,factor in [("fraction",1),("percent",100)]}
closest = min(errors, key=errors.get)
state = {"source": "quantitative MRI volume columns and independently extracted gene expression rows",
         "patients": paired.patient.nunique(), "records": len(paired),
         "max_endpoint_difference": float((paired.y-paired.reference_y).abs().max()),
         "max_clinical_input_difference": float((paired.C_cyt_late-paired.score).abs().max()),
         "reported_change_units": closest, "reported_change_maximum_discrepancy": errors[closest],
         "reserved_outcomes_used": False}
(Path(study["report_directory"]) / "endpoint_verification.json").write_text(json.dumps(state, indent=2) + "\n")
print(json.dumps(state))
