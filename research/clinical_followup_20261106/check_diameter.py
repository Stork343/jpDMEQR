import json
from pathlib import Path

import numpy as np
import pandas as pd
from openpyxl import load_workbook
from scipy.stats import rankdata

from panel import study_config

study = study_config("ispy_diameter")
inputs = Path(study["input_directory"])
src = Path(study["source_directory"])
roles = pd.read_csv(inputs / "roles.csv", dtype={"source_id": str}).query("role == 'discovery'").set_index("source_id")
matrix = pd.read_csv(src / "ispy1_panel_expression.csv").set_index("gene")
program = ["MKI67", "BIRC5", "CCNB1"]
book = load_workbook(src / "ispy1_clinical.xlsx", read_only=True, data_only=True)
values = book["TCIA Patient Clinical Subset"].values
header = list(next(values))
endpoints = ["MRI LD Baseline", "MRI LD 1-3dAC", "MRI LD InterReg", "MRI LD PreSurg"]
reconstructed = []
for record in values:
    origin = str(record[header.index("SUBJECTID")])
    if origin not in roles.index:
        continue
    baseline, *follow = [record[header.index(field)] for field in endpoints]
    if baseline is None or float(baseline) <= 0:
        continue
    x = matrix[origin]
    observed = x.notna()
    rank = pd.Series(rankdata(x[observed], method="average") / observed.sum(), index=x.index[observed])
    score = rank.reindex(program).mean(skipna=False)*10 if observed.mean() >= .8 else np.nan
    for period, value in enumerate(follow, 1):
        if value is not None and float(value) >= 0:
            key = "MRI_"+roles.at[origin,"patient"]+"_"+str(period)
            reconstructed.append((key, np.log1p(float(value)/float(baseline)), score*(period-1)/2))
book.close()
reference = pd.DataFrame(reconstructed, columns=["sample_id","ref_y","ref_score"])
rows = pd.read_csv(inputs/"discovery.csv")
features = pd.read_csv(inputs/"features.csv.gz", usecols=["sample_id","C_proliferation_late"])
both = reference.merge(rows, on="sample_id", validate="one_to_one").merge(features, on="sample_id", validate="one_to_one")
assert len(both) == len(reference) == len(rows)
np.testing.assert_allclose(both.y, both.ref_y, atol=1e-12, rtol=0)
np.testing.assert_allclose(both.C_proliferation_late, both.ref_score, atol=1e-12, rtol=0, equal_nan=True)
state = {"patients": both.strain_id.nunique(), "records": len(both),
         "source": "TCIA clinical workbook and independent within-array ranking of the common GEO measurements",
         "maximum_endpoint_difference": float(abs(both.y-both.ref_y).max()),
         "maximum_contrast_difference": float(abs(both.C_proliferation_late-both.ref_score).max()),
         "reserved_outcomes_used": False}
(Path(study["report_directory"])/"endpoint_verification.json").write_text(json.dumps(state, indent=2)+"\n")
print(json.dumps(state))
