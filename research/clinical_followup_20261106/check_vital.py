from pathlib import Path
import json

import duckdb
import numpy as np
import pandas as pd

from panel import study_config

study = study_config("vitaldb_crp")
inputs = Path(study["input_directory"])
source = Path(study["source_directory"])
roles = pd.read_csv(inputs / "roles.csv")
assert roles.patient.is_unique and roles.caseid.is_unique
db = duckdb.connect()
db.register("roles", roles)
db.read_csv(str(source / "cases.csv.gz"), all_varchar=True).create_view("cases")
db.read_csv(str(source / "labs.csv.gz"), all_varchar=True).create_view("labs")
index = db.sql("""
    select r.patient,r.caseid,greatest(cast(c.caseend as double),cast(c.aneend as double),cast(c.opend as double)) as finish,
      cast(c.dis as double) as discharge,case when c.approach='Open' then (cast(c.opend as double)-cast(c.opstart as double))/3600-3
        else 0 end as contrast
      from roles r join cases c on r.caseid=cast(c.caseid as integer)
      where r.role='discovery' and cast(c.subjectid as integer)=r.patient
""")
index.create_view("cohort")
outcome = db.sql("""
    with records as (select c.patient,c.caseid,c.contrast,cast(l.dt as double) as t,(t-c.finish)/3600 as hours,
         median(cast(l.result as double)) as value from cohort c join labs l on c.caseid=cast(l.caseid as integer)
         where l.name='crp' and cast(l.result as double)>=0 and t<=c.discharge and hours>12 and hours<=84
         group by c.patient,c.caseid,c.contrast,t,c.finish),
    windows as (select *, case when hours<=36 then 1 when hours<=60 then 2 else 3 end as period from records)
    select * from windows qualify row_number() over(partition by patient,period order by abs(hours-period*24),t)=1
""").df()
outcome["sample_id"] = "VITAL_" + outcome.patient.astype(str) + "_" + outcome.period.astype(str)
rows = pd.read_csv(inputs / "discovery.csv")
x = pd.read_csv(inputs / "features.csv.gz", usecols=["sample_id", "C_open_duration"])
paired = rows.merge(outcome, on="sample_id", validate="one_to_one").merge(x, on="sample_id", validate="one_to_one")
assert len(paired) == len(rows) == len(outcome)
np.testing.assert_allclose(paired.y, paired.value, atol=1e-10)
np.testing.assert_allclose(paired.C_open_duration, paired.contrast, atol=1e-10)
state = {"source": "independent DuckDB reconstruction of the historical VitalDB case and laboratory tables",
         "patients": outcome.patient.nunique(), "records": len(outcome),
         "maximum_endpoint_difference": float(abs(paired.y-paired.value).max()),
         "maximum_clinical_input_difference": float(abs(paired.C_open_duration-paired.contrast).max()),
         "reserved_outcomes_used": False}
(Path(study["report_directory"]) / "endpoint_verification.json").write_text(json.dumps(state, indent=2) + "\n")
print(json.dumps(state))
