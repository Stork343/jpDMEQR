from pathlib import Path
from zipfile import ZipFile
import io
import json

import duckdb
import numpy as np
import pandas as pd

from panel import study_config

study = study_config("labels_t1dx_monitoring")
inputs = Path(study["input_directory"])
source = Path(study["source_directory"])
roles = pd.read_csv(inputs/"roles.csv",dtype={"patient":str})
assert roles.patient.is_unique and roles.groupby("SiteID").role.nunique().eq(1).all()
db = duckdb.connect()
db.register("roles",roles)
with ZipFile(source/"registry.zip") as archive:
    for view,suffix in [("visits","/VisitsLong.txt"),("lab","/Hba1cLong.txt")]:
        path, = [f for f in archive.namelist() if f.endswith(suffix)]
        db.read_csv(io.BytesIO(archive.read(path)),delimiter="|",all_varchar=True).to_table(view)
verified = db.sql("""
    with cohort as (
      select r.patient,greatest(0,r.ExamDtMnC) as landmark,
             (cast(v.BGTestAvgNumPtRep as double)-cast(v.BGTestAvgNumMeter as double))/2 as gap,
             case v.Pt_Gender when 'F' then 1 when 'M' then 0 when 'T' then 0 end as female
      from roles r join visits v on r.patient=v.PtID and r.Visit=v.Visit and r.ExamDtMnC=cast(v.ExamDtMnC as double)
      where r.role='discovery'), assays as (
      select c.patient,cast(l.HbA1cImputeDtMnC as double) as t,t-c.landmark as delay,
             median(cast(l.HbA1c as double)) as y,c.gap,c.female from lab l join cohort c on l.PtID=c.patient
      where l.HbA1c is not null and l.HbA1cImputeDtMnC is not null group by c.patient,t,c.landmark,c.gap,c.female),
    eligible as (select distinct patient from assays where delay between -6 and 0),
    windows as (select a.*,case when delay<=9 then 1 when delay<=15 then 2 else 3 end as period
      from assays a join eligible using(patient) where delay>3 and delay<=27)
    select * from windows qualify row_number() over(partition by patient,period order by
      abs(delay-(case period when 1 then 6 when 2 then 12 else 24 end)),t)=1
""").df()
verified["sample_id"] = "DX_"+verified.patient+"_"+verified.period.astype(str)
rows = pd.read_csv(inputs/"discovery.csv")
features = pd.read_csv(inputs/"features.csv.gz",usecols=["sample_id","C_check_gap","C_female"])
paired = rows.merge(verified,on="sample_id",validate="one_to_one",suffixes=("_stored","_sql")).merge(features,on="sample_id",validate="one_to_one")
assert len(paired)==len(rows)==len(verified)
np.testing.assert_allclose(paired.y_stored,paired.y_sql,atol=1e-12)
np.testing.assert_allclose(paired.C_check_gap,paired.gap,atol=1e-12)
np.testing.assert_allclose(paired.C_female,paired.female,atol=1e-12,equal_nan=True)
state = {"source":"independent DuckDB reconstruction of the public longitudinal visit and HbA1c tables",
         "patients":paired.patient.nunique(),"records":len(paired),"maximum_endpoint_difference":float(abs(paired.y_stored-paired.y_sql).max()),
         "maximum_monitoring_input_difference":float(abs(paired.C_check_gap-paired.gap).max()),"sex_codes_verified":True,
         "reserved_outcomes_used":False}
(Path(study["report_directory"])/"endpoint_verification.json").write_text(json.dumps(state,indent=2)+"\n")
print(json.dumps(state))
