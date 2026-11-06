import json
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

from panel import study_config

study = study_config("gaitrec_propulsion")
inputs = Path(study["input_directory"])
source = Path(study["source_directory"])
audit = pd.read_csv(inputs / "time_audit.csv")
rows = pd.read_csv(inputs / "discovery.csv")
features = pd.read_csv(inputs / "features.csv.gz", usecols=["sample_id", "C_propulsion_delay"])
roles = pd.read_csv(inputs / "roles.csv").set_index("patient")
assert roles.loc[audit.patient, "role"].eq("discovery").all()
database = duckdb.connect()
database.register("audit", audit)
for side in ["left", "right"]:
    path = source / f"GRF_F_AP_PRO_{side}.csv"
    database.read_csv(str(path), all_varchar=True).create_view("source_" + side)
all_columns = [f"F_AP_PRO_{i}" for i in range(1, 102)]
nonnumeric = " or ".join(f"try_cast({c} as double) is null" for c in all_columns)
peak = "greatest(0," + ",".join(f"cast(F_AP_PRO_{t + 1} as double)" for t in range(50, 96)) + ")"
mass = " + ".join(f"greatest(0,cast(F_AP_PRO_{t+1} as double))" for t in range(30, 96))
moment = " + ".join(f"{t}*greatest(0,cast(F_AP_PRO_{t+1} as double))" for t in range(30, 96))
for side in ["left", "right"]:
    database.sql(f"""
        select cast(SUBJECT_ID as bigint) as patient,cast(SESSION_ID as bigint) as session_id,TRIAL_ID,
               {peak} as peak,({moment})/nullif(({mass}),0) as centroid
          from source_{side} s where exists (select 1 from audit a where a.patient=cast(s.SUBJECT_ID as bigint)
                 and cast(s.SESSION_ID as bigint) in (a.baseline_session,a.follow_session)) and not ({nonnumeric})
    """).create_view(side + "_ap")
independent = database.sql("""
    with paired as (select l.patient,l.session_id,l.TRIAL_ID,l.peak as left_peak,r.peak as right_peak,
                    l.centroid-r.centroid as timing from left_ap l join right_ap r using(patient,session_id,TRIAL_ID))
    select a.sample_id,100*median(case when a.affected='left' then p.left_peak else p.right_peak end) as y,
           count(*) as paired_trials from audit a join paired p on p.patient=a.patient and p.session_id=a.follow_session
           group by a.sample_id
""").df()
paired = rows.merge(independent, on="sample_id", validate="one_to_one", suffixes=("_saved", "_sql"))
assert len(paired) == len(rows) == len(audit)
error = float(abs(paired.y_sql - paired.y_saved).max())
assert error < 1e-9
timing = database.sql("""
    select a.sample_id,median((case when a.affected='left' then 1 else -1 end)*(l.centroid-r.centroid))/10 as value
        from audit a join left_ap l on l.patient=a.patient and l.session_id=a.baseline_session
        join right_ap r on r.patient=l.patient and r.session_id=l.session_id and r.TRIAL_ID=l.TRIAL_ID
        group by a.sample_id
""").df()
timing = features.merge(timing, on="sample_id", validate="one_to_one")
np.testing.assert_allclose(timing.C_propulsion_delay, timing.value, atol=1e-9, equal_nan=True)
state = {"source": "independent SQL aggregation of published bilateral AP signals", "records": len(paired),
         "patients": audit.patient.nunique(), "maximum_outcome_difference": error,
         "maximum_timing_difference": float(abs(timing.C_propulsion_delay - timing.value).max()),
         "reserved_outcomes_used": False}
(Path(study["report_directory"]) / "endpoint_verification.json").write_text(json.dumps(state, indent=2) + "\n")
print(json.dumps(state))
