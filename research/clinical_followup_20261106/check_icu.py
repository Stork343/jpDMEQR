from pathlib import Path
import json

import duckdb
import numpy as np
import pandas as pd

from panel import study_config

for name, lab, unit, target in [("eicu_lactate", "lactate", "mmol/L", "dsi"),
                                ("eicu_creatinine", "creatinine", "mg/dL", "pulse")]:
    study = study_config(name)
    inputs = Path(study["input_directory"])
    source = Path(study["source_directory"])
    roles = pd.read_csv(inputs / "roles.csv", dtype={"patient": str})
    assert roles.patient.is_unique and roles.groupby("hospitalid").role.nunique().eq(1).all()
    database = duckdb.connect()
    database.register("roles", roles)
    for table in ["lab", "patient", "vitalPeriodic", "vitalAperiodic"]:
        database.read_csv(str(source / f"{table}.csv.gz"), all_varchar=True, header=True).create_view(table)
    observed = database.sql(f"""
        with samples as (
          select r.patient, r.patientunitstayid as stay, cast(l.labresultoffset as integer) as t,
            cast(l.labresultrevisedoffset as integer) as entered, cast(l.labresult as double) as val
          from lab l join roles r on cast(l.patientunitstayid as integer)=r.patientunitstayid
          join patient p on p.patientunitstayid=l.patientunitstayid
          where r.role='discovery' and l.labname='{lab}' and l.labmeasurenamesystem='{unit}'
            and not regexp_matches(coalesce(l.labresulttext,''),'[<>≤≥]') and val>0
            and t between 0 and least(2880, cast(p.unitdischargeoffset as integer))
        ), qualified as (
          select stay from samples where t<=720 and entered<=720 group by stay having max(val)>={2 if lab == 'lactate' else 0}
        ), future as (
          select patient,stay,t,((t-721)//720)+1 as period,median(val) as y from samples
          join qualified using(stay) where t>720 group by patient,stay,t
        )
        select * from future qualify row_number() over
          (partition by stay,period order by abs(t-(period*720+360)),t)=1
    """).df()
    observed["sample_id"] = name + "_" + observed.patient + "_" + observed.period.astype(str)
    rows = pd.read_csv(inputs / "discovery.csv")
    paired = observed.merge(rows, on="sample_id", suffixes=("_sql", "_saved"), validate="one_to_one")
    assert len(paired) == len(rows) == len(observed)
    maximum = float(abs(paired.y_sql - paired.y_saved).max())
    assert maximum < 1e-9
    vital = database.sql("""
        with invasive as (
          select r.patientunitstayid as stay, greatest(0,(cast(observationoffset as integer)-1)//60) as hour,
            median(try_cast(heartrate as double)) filter (where try_cast(heartrate as double)>0 and try_cast(heartrate as double)<=300) as hr,
            median(try_cast(systemicsystolic as double)) filter (where try_cast(systemicsystolic as double)>0 and try_cast(systemicsystolic as double)<=350) as s,
            median(try_cast(systemicdiastolic as double)) filter (where try_cast(systemicdiastolic as double)>0 and try_cast(systemicdiastolic as double)<=300) as d
          from vitalPeriodic v join roles r on cast(v.patientunitstayid as integer)=r.patientunitstayid
          join patient p on p.patientunitstayid=v.patientunitstayid
          where r.role='discovery' and cast(observationoffset as integer) between 0 and least(720,cast(p.unitdischargeoffset as integer))
          group by stay,hour
        ), cuff as (
          select r.patientunitstayid as stay, greatest(0,(cast(observationoffset as integer)-1)//60) as hour,
            median(try_cast(noninvasivesystolic as double)) filter (where try_cast(noninvasivesystolic as double)>0 and try_cast(noninvasivesystolic as double)<=350) as s,
            median(try_cast(noninvasivediastolic as double)) filter (where try_cast(noninvasivediastolic as double)>0 and try_cast(noninvasivediastolic as double)<=300) as d
          from vitalAperiodic v join roles r on cast(v.patientunitstayid as integer)=r.patientunitstayid
          join patient p on p.patientunitstayid=v.patientunitstayid
          where r.role='discovery' and cast(observationoffset as integer) between 0 and least(720,cast(p.unitdischargeoffset as integer))
          group by stay,hour
        ), measures as (
          select coalesce(i.stay,c.stay) as stay,i.hr,
            case when i.s>=i.d then i.s when c.s>=c.d then c.s end as systolic,
            case when i.s>=i.d then i.d when c.s>=c.d then c.d end as diastolic
          from invasive i full join cuff c on i.stay=c.stay and i.hour=c.hour
        )
        select stay, case when count(hr/diastolic)>=2 then median(hr/diastolic) end as dsi,
          case when count(systolic-diastolic)>=2 then median((systolic-diastolic)/10) end as pulse
          from measures group by stay
    """).df()
    inputs_x = pd.read_csv(inputs / "features.csv.gz", usecols=["sample_id", study["inference_target"]])
    check = paired.merge(vital, on="stay", how="left", validate="many_to_one").merge(inputs_x, on="sample_id", validate="one_to_one")
    np.testing.assert_allclose(check[target], check[study["inference_target"]], atol=1e-10, equal_nan=True)
    state = {"source": "independent DuckDB reconstruction from published compressed tables",
             "patients": observed.patient.nunique(), "records": len(observed), "maximum_endpoint_difference": maximum,
             "maximum_clinical_input_difference": float(abs(check[target] - check[study["inference_target"]]).max()),
             "reserved_outcomes_used": False}
    (Path(study["report_directory"]) / "endpoint_verification.json").write_text(json.dumps(state, indent=2) + "\n")
    print(json.dumps(state), flush=True)
    database.close()
