from csv import DictReader
from io import TextIOWrapper
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd

from panel import study_config

study = study_config("physionet_ph")
inputs = Path(study["input_directory"])
src = Path(study["source_directory"])
cohort = pd.read_csv(src/"catalog.csv").query("site == 'A'")
roles = pd.read_csv(inputs/"roles.csv")
assert set(cohort.patient) == set(roles.query("role == 'discovery'").patient)
series = []
for batch, group in cohort.groupby("batch"):
    with ZipFile(src/f"A_{batch:03}.zip") as archive:
        for patient in group.patient:
            with archive.open(patient+".psv") as stream:
                visits = [(int(v["ICULOS"]),float(v["Age"]),float(v["pH"]),float(v["PaCO2"]),float(v["Resp"]))
                          for v in DictReader(TextIOWrapper(stream), delimiter="|") if int(v["ICULOS"]) <= 48]
            baseline = [v for v in visits if 1 <= v[0] <= 12]
            ph = [v for v in baseline if 6.5 <= v[2] <= 8]
            if not ph or not any(v[1]>=18 for v in baseline):
                continue
            pairs = [v for v in baseline if v[3]>0 and v[4]>0]
            pair = max(pairs, key=lambda v: v[0]) if pairs else None
            contrast = (pair[3]-40)*(pair[4]-20)/50 if pair else np.nan
            for period, center in enumerate([18,30,42], 1):
                valid = [v for v in visits if center-5 <= v[0] <= center+6 and 6.5 <= v[2] <= 8]
                if valid:
                    observed = min(valid, key=lambda v: (abs(v[0]-center),v[0]))
                    series.append((patient+"_pH_"+str(period),observed[2],contrast,max(v[0] for v in ph),
                                   observed[0],pair[0] if pair else np.nan))
reference = pd.DataFrame(series, columns=["sample_id","ref_y","ref_c","base","hour","pair"])
rows = pd.read_csv(inputs/"discovery.csv")
features = pd.read_csv(inputs/"features.csv.gz", usecols=["sample_id","C_ventilation"])
audit = pd.read_csv(inputs/"time_audit.csv")
both = reference.merge(rows,on="sample_id",validate="one_to_one").merge(features,on="sample_id",validate="one_to_one")
both = both.merge(audit, on="sample_id",validate="one_to_one")
assert len(both) == len(rows) == len(reference) == len(audit)
for original, rebuilt in [(both.y,both.ref_y),(both.C_ventilation,both.ref_c),(both.baseline_hour,both.base),
                          (both.ICULOS,both.hour),(both.pair_hour,both.pair)]:
    np.testing.assert_allclose(original, rebuilt, atol=1e-12, rtol=0, equal_nan=True)
state = {"source": "second reading of the original hospital A PSV files using a separate patient and visit selector",
         "patients": both.strain_id.nunique(), "records": len(both),
         "maximum_endpoint_difference": float(abs(both.y-both.ref_y).max()),
         "maximum_contrast_difference": float(abs(both.C_ventilation-both.ref_c).max()),
         "hours_verified": True, "reserved_outcomes_used": False}
(Path(study["report_directory"])/"endpoint_verification.json").write_text(json.dumps(state, indent=2)+"\n")
print(json.dumps(state))
