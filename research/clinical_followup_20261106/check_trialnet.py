from pathlib import Path
import json
import re

import duckdb
import numpy as np
import pandas as pd

from panel import study_config

study = study_config("trialnet_cpeptide")
src = Path(study["source_directory"])
inputs = Path(study["input_directory"])
roles = pd.read_csv(inputs/"roles.csv",dtype={"patient":str})
cohort = roles.query("role=='discovery'").set_index("patient")
assert not cohort.study.isin(["TN02","START"]).any()
factor = 1000/3020.29/120
floor = 6*factor
minutes = [0,15,30,60,90,120]
db = duckdb.connect()
observations, counts = [], []
for tag,prefix in [("controls","")]:
    members = cohort.index[cohort.source.eq(tag)]
    f = pd.read_csv(src/(prefix+"cpeptide.csv"),dtype={"participant_id":str},
                    usecols=["participant_id","cpeptide_study_day",*["pep"+str(t) for t in minutes]])
    f = f.loc[f.participant_id.isin(members)].copy()
    curve = f[["pep"+str(t) for t in minutes]]
    valid = curve.notna().all(axis=1) & curve.ge(0).all(axis=1)
    f["amount"] = np.trapz(curve.clip(lower=.05).to_numpy(),minutes,axis=1)*factor
    f.loc[~valid,"amount"] = np.nan
    observations.append(f[["participant_id","cpeptide_study_day","amount"]].rename(columns={
      "participant_id":"patient","cpeptide_study_day":"day"}))
    b = pd.read_csv(src/(prefix+"cbc.csv"),dtype={"participant_id":str},
                    usecols=["participant_id","cbc_study_day","neutrophils_abs","lymphocytes_abs"])
    counts.append(b.loc[b.participant_id.isin(members)].rename(columns={"participant_id":"patient","cbc_study_day":"day"}))

reference = pd.read_csv(src/"cpeptide.csv",dtype={"participant_id":str},
  usecols=["study","participant_id","cpeptide_visit_number","auc2hr_ng_ml","auc2hr_nmol_l_min",
           *["pep"+str(t) for t in minutes]])
reference = reference.loc[reference.study.eq("TN09")].copy()
reference["weeks"] = reference.cpeptide_visit_number.map({-1:0,4:12,8:28,14:52,21:80,27:104,28:130,29:156})
published = pd.read_csv(src/"tn09_cpeptide.txt",sep="\t",dtype={"id":str},
                        usecols=["id","weeks","aucResults",*["X"+str(t) for t in minutes]])
published["patient"] = "TN09_"+published.id
curves = published[["X"+str(t) for t in minutes]].apply(pd.to_numeric,errors="raise")
available = curves.notna().all(axis=1) & curves.ge(0).all(axis=1)
published["report_integral"] = np.trapz(curves.clip(lower=.05).to_numpy(),minutes,axis=1)
published.loc[~available,"report_integral"] = np.nan
unit = reference.merge(published,left_on=["participant_id","weeks"],right_on=["patient","weeks"],validate="one_to_one")
assert len(unit)>200
curves = unit[["pep"+str(t) for t in minutes]]
complete = curves.notna().all(axis=1) & curves.ge(0).all(axis=1)
integral = np.trapz(curves.loc[complete].clip(lower=.05).to_numpy(),minutes,axis=1)
np.testing.assert_allclose(integral,unit.loc[complete,"report_integral"],rtol=0,atol=1e-8)
np.testing.assert_allclose(unit.auc2hr_nmol_l_min,unit.auc2hr_ng_ml*factor,rtol=0,atol=1e-12)
members = cohort.index[cohort.source.eq("abatacept")]
b = pd.read_csv(src/"tn09_cbc.csv",dtype={"participant_id":str},
                usecols=["participant_id","cbc_study_day","neutrophils_abs","lymphocytes_abs","day_0"])
b = b.loc[b.participant_id.isin(members)].rename(columns={"participant_id":"patient","cbc_study_day":"day"})
counts.append(b.drop(columns="day_0"))
origin = b[["patient","day_0"]].drop_duplicates().set_index("patient").day_0
clinical = pd.read_csv(src/"tn09_clinical.txt",sep="\t",dtype={"Random_ID":str},na_values=["."])
clinical["patient"] = "TN09_"+clinical.Random_ID
clinical = clinical.loc[clinical.patient.isin(members)]
columns = list(clinical.columns)
week_of_visit = {"Screening":0,"Day 84":12,"Day 196":28,"Day 364":52,"Day 560":80,"Day 728":104,
                 "2 years, 6 months":130,"3 years":156}
draws = []
for field in [c for c in columns if re.fullmatch(r"AUC(?:\.\d+)?",c)]:
    index = columns.index(field)
    visit = max(i for i in range(index) if re.fullmatch(r"Visit(?:\.\d+)?",columns[i]))
    sample_date = max(i for i in range(visit,index) if re.fullmatch(r"Draw Date(?:\.\d+)?",columns[i]))
    curve_names = [next(columns[i] for i in range(visit,index)
                        if re.fullmatch(rf"PEP{t}(?:\.\d+)?",columns[i])) for t in minutes]
    curve = clinical[curve_names].replace({"<.05":.05,"<0.05":.05}).apply(pd.to_numeric,errors="raise")
    valid = curve.notna().all(axis=1) & curve.ge(0).all(axis=1)
    amount = pd.Series(np.trapz(curve.clip(lower=.05).to_numpy(),minutes,axis=1),index=clinical.index).where(valid)
    frame = pd.DataFrame({"patient":clinical.patient,"weeks":clinical.iloc[:,visit].map(week_of_visit),
                          "draw":clinical.iloc[:,sample_date],"enrollment":clinical.patient.map(origin),
                          "wide_integral":amount})
    draws.append(frame)
draws = pd.concat(draws,ignore_index=True).merge(published[["patient","weeks","report_integral"]],
                                              on=["patient","weeks"],validate="one_to_one")
complete = draws.wide_integral.notna() & draws.report_integral.notna()
np.testing.assert_allclose(draws.loc[complete,"wide_integral"],draws.loc[complete,"report_integral"],rtol=0,atol=1e-8)
draws = draws.loc[complete]
db.register("draws",draws.dropna(subset=["draw","enrollment","report_integral"]))
active = db.sql("""
  select patient,date_diff('day',strptime(enrollment,['%Y-%m-%d','%m/%d/%y','%m/%d/%Y']),
      strptime(draw,['%Y-%m-%d','%m/%d/%y','%m/%d/%Y'])) as day,report_integral as integral from draws
""").df()
active["amount"] = active.integral*factor
observations.append(active[["patient","day","amount"]])
measurements = pd.concat(observations,ignore_index=True).dropna()
measurements["amount"] = np.maximum(measurements.amount,floor)
db.register("measurements",measurements)
expected = db.sql("""
  with assays as (select patient,day,median(amount) as amount from measurements group by patient,day),
  baseline as (select patient,day as baseline_day,amount as baseline_amount from assays
    where day between -60 and 30 qualify row_number() over(partition by patient order by day desc)=1),
  visits as (select a.*,b.baseline_day,b.baseline_amount,
      case when day between 120 and 240 then 1 when day between 300 and 450 then 2
           when day between 600 and 840 then 3 end as period
    from assays a join baseline b using(patient))
  select *,ln(1+amount) as reference_y from visits where period is not null
  qualify row_number() over(partition by patient,period order by
      abs(day-case period when 1 then 182 when 2 then 365 else 730 end),day)=1
""").df()
blood = pd.concat(counts,ignore_index=True).dropna(subset=["day"])
db.register("blood",blood)
immune = db.sql("""
  with assays as (select patient,day,median(neutrophils_abs) as neutrophils,median(lymphocytes_abs) as lymphocytes
    from blood where day between -30 and 30 group by patient,day)
  select patient,day as cbc_day,ln(neutrophils/lymphocytes) as log_ratio from assays
    where neutrophils>0 and lymphocytes>0 qualify row_number() over(partition by patient order by day desc)=1
""").df()
expected = expected.merge(immune,on="patient",how="left",validate="many_to_one")
expected["sample_id"] = expected.patient+"_CP_"+expected.period.astype(str)
rows = pd.read_csv(inputs/"discovery.csv")
features = pd.read_csv(inputs/"features.csv.gz",usecols=["sample_id","C_nlr_late","C_baseline"])
both = expected.merge(rows,on="sample_id",validate="one_to_one").merge(features,on="sample_id",validate="one_to_one")
audit = pd.read_csv(inputs/"time_audit.csv")
both = both.merge(audit[["sample_id","day","baseline_day","cbc_day"]],on="sample_id",suffixes=("_new","_saved"),validate="one_to_one")
assert len(both)==len(rows)==len(expected)==len(audit),(len(both),len(rows),len(expected),len(audit))
assert set(both.patient) <= set(cohort.index)
for column in ["day","baseline_day","cbc_day"]:
    np.testing.assert_allclose(both[column+"_new"],both[column+"_saved"],atol=1e-12,rtol=0,equal_nan=True)
np.testing.assert_allclose(both.y,both.reference_y,atol=1e-12,rtol=0)
np.testing.assert_allclose(both.C_baseline,np.log1p(both.baseline_amount),atol=1e-12,rtol=0)
target = both.log_ratio*(both.period-1)/2
np.testing.assert_allclose(both.C_nlr_late,target,atol=1e-12,rtol=0,equal_nan=True)
state = {"source":"independent SQL visit selection and date arithmetic, author-normalized and integral C-peptide tables",
         "patients":both.patient.nunique(),"records":len(both),"overlapping_unit_checks":len(unit),
         "maximum_endpoint_difference":float(abs(both.y-both.reference_y).max()),
         "maximum_contrast_difference":float(abs(both.C_nlr_late-target).max()),
         "reserved_outcomes_used":False}
(Path(study["report_directory"])/"endpoint_verification.json").write_text(json.dumps(state,indent=2)+"\n")
print(json.dumps(state))
