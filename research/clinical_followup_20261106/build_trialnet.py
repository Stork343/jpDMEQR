from pathlib import Path
import re

import numpy as np
import pandas as pd

from panel import finish_panel, study_config

study = study_config("trialnet_cpeptide")
src = Path(study["source_directory"])
inputs = Path(study["input_directory"])
roles = pd.read_csv(inputs/"roles.csv",dtype={"patient":str}).query("role=='discovery'").set_index("patient")
assert not roles.study.isin(["TN02","START"]).any()
conversion = 1000/3020.29/120
floor = .05*1000/3020.29
minutes = np.array([0,15,30,60,90,120])
integration_weights = np.array([7.5,15,22.5,30,30,15])


def reportable_measurements(frame):
    output = pd.DataFrame(index=frame.index)
    for column in frame:
        text = frame[column].astype("string").str.strip()
        bounded = text.str.fullmatch(r"<\s*(?:\d+(?:\.\d*)?|\.\d+)",na=False)
        limits = pd.to_numeric(text.where(bounded).str.replace("<","",regex=False),errors="raise")
        exact = pd.to_numeric(text.mask(bounded),errors="raise")
        exact = exact.where(exact.ge(0))
        output[column] = exact.clip(lower=.05).mask(bounded & limits.le(.05),.05).astype(float)
    return output


parts, blood = [], []
for tag, prefix in [("controls","")]:
    ids = roles.index[roles.source.eq(tag)]
    frame = pd.read_csv(src/(prefix+"cpeptide.csv"),dtype={"participant_id":str},
                        usecols=["participant_id","cpeptide_study_day","auc2hr_ng_ml","auc2hr_nmol_l_min",
                                 *["pep"+str(t) for t in minutes]])
    frame = frame.loc[frame.participant_id.isin(ids)].copy()
    np.testing.assert_allclose(frame.auc2hr_ng_ml*conversion,frame.auc2hr_nmol_l_min,rtol=0,atol=1e-12,equal_nan=True)
    curve = reportable_measurements(frame[["pep"+str(t) for t in minutes]])
    valid = curve.notna().all(axis=1) & curve.ge(0).all(axis=1)
    frame["integral"] = curve.clip(lower=.05).mul(integration_weights,axis=1).sum(axis=1).where(valid)
    frame = frame.rename(columns={"participant_id":"patient","cpeptide_study_day":"day"})
    parts.append(frame[["patient","day","integral"]])
    fields = ["participant_id","cbc_study_day","neutrophils_abs","lymphocytes_abs","wbc"]
    counts = pd.read_csv(src/(prefix+"cbc.csv"),usecols=fields,dtype={"participant_id":str})
    blood.append(counts.loc[counts.participant_id.isin(ids)].rename(columns={"participant_id":"patient","cbc_study_day":"day"}))

cbc = pd.read_csv(src/"tn09_cbc.csv",dtype={"participant_id":str},
                  usecols=["participant_id","cbc_study_day","neutrophils_abs","lymphocytes_abs","wbc","day_0"])
cbc = cbc.rename(columns={"participant_id":"patient","cbc_study_day":"day"})
ids = roles.index[roles.source.eq("abatacept")]
assert cbc.loc[cbc.patient.isin(ids)].groupby("patient").day_0.nunique().le(1).all()
anchors = cbc.loc[cbc.patient.isin(ids)].drop_duplicates("patient").set_index("patient").day_0
blood.append(cbc.loc[cbc.patient.isin(ids)].drop(columns="day_0"))
clinical = pd.read_csv(src/"tn09_clinical.txt",sep="\t",dtype={"Random_ID":str},na_values=["."])
clinical["patient"] = "TN09_"+clinical.Random_ID
clinical = clinical.loc[clinical.patient.isin(ids)].set_index("patient")
headers = list(clinical.columns)
for column in [c for c in headers if re.fullmatch(r"AUC(?:\.\d+)?",c)]:
    position = headers.index(column)
    visit_position = max(i for i in range(position) if re.fullmatch(r"Visit(?:\.\d+)?",headers[i]))
    draw_position = max(i for i in range(visit_position,position)
                        if re.fullmatch(r"Draw Date(?:\.\d+)?",headers[i]))
    date_column = headers[draw_position]
    day = (pd.to_datetime(clinical[date_column],format="mixed")-pd.to_datetime(anchors,format="mixed")).dt.days
    curve_columns = [next(headers[i] for i in range(visit_position,position)
                         if re.fullmatch(rf"PEP{t}(?:\.\d+)?",headers[i])) for t in minutes]
    curve = reportable_measurements(clinical[curve_columns])
    valid = curve.notna().all(axis=1) & curve.ge(0).all(axis=1)
    amount = curve.clip(lower=.05).mul(integration_weights,axis=1).sum(axis=1).where(valid)
    parts.append(pd.DataFrame({"day":day,"integral":amount}).rename_axis("patient").reset_index())
assays = pd.concat(parts,ignore_index=True).dropna(subset=["day","integral"])
assert set(assays.patient) <= set(roles.index) and np.isfinite(assays.integral).all()
assert assays.integral.ge(0).all()
assays = assays.groupby(["patient","day"],as_index=False).integral.median()
assays["mean_cpeptide"] = assays.integral*conversion
baseline = assays.loc[assays.day.between(-60,30)].sort_values(["patient","day"]).drop_duplicates("patient",keep="last").set_index("patient")
future = assays.loc[assays.patient.isin(baseline.index)].copy()
future["period"] = np.select([future.day.between(120,240),future.day.between(300,450),future.day.between(600,840)],[1,2,3],0)
future = future.loc[future.period.gt(0)]
future["distance"] = abs(future.day-future.period.map({1:182,2:365,3:730}))
future = future.sort_values(["patient","period","distance","day"]).drop_duplicates(["patient","period"])
assert not future.empty
patients = pd.Index(sorted(future.patient.unique()))
counts = pd.concat(blood,ignore_index=True).dropna(subset=["day"])
assert set(counts.patient) <= set(roles.index)
counts = counts.loc[counts.day.between(-30,30)].groupby(["patient","day"],as_index=False).median(numeric_only=True)
paired = counts.loc[counts.neutrophils_abs.gt(0) & counts.lymphocytes_abs.gt(0)].sort_values("day").drop_duplicates("patient",keep="last").set_index("patient")
last = counts.sort_values("day").drop_duplicates("patient",keep="last").set_index("patient")
meta = roles.loc[patients]
features = pd.DataFrame(index=patients)
features["C_log_age"] = np.log(meta.age_years)
sex = meta.sex.str.lower().map({"female":1,"male":0,"f":1,"m":0})
assert sex.notna().all()
features["C_female"] = sex
features["C_baseline"] = np.log1p(baseline.mean_cpeptide)
features["C_days_since_test"] = 30-baseline.day
features["C_nlr"] = np.log(paired.neutrophils_abs/paired.lymphocytes_abs)
for variable in ["neutrophils_abs","lymphocytes_abs","wbc"]:
    features["C_"+variable] = np.log1p(last[variable].where(last[variable].ge(0)))
for therapy in ["rituximab","abatacept"]:
    features["C_"+therapy] = meta.treatment.eq(therapy).astype(int)
hba = pd.read_csv(src/"patients.csv",usecols=["participant_id","baseline_hba1c"],dtype={"participant_id":str}).set_index("participant_id").baseline_hba1c
hba = hba.reindex(roles.index)
hba.loc[clinical.index] = pd.to_numeric(clinical["HbA1c"],errors="raise")
features["C_hba1c"] = hba

rna = pd.read_csv(src/"baseline_rna.csv",dtype={"patient":str})
rna = rna.loc[rna.patient.isin(patients)].copy()
rna["distance"] = abs(rna.day)
day = rna.sort_values(["patient","distance","day"]).drop_duplicates("patient").set_index("patient").day
rna = rna.loc[rna.day.eq(rna.patient.map(day))]
tables = {"controls":("counts.csv.gz",","),"rituximab":("tn05_counts.gz","\t"),"abatacept":("tn09_counts.gz",",")}
genes, matrices = None, {}
for tag,(filename,separator) in tables.items():
    columns = pd.read_csv(src/filename,sep=separator,nrows=0).columns
    gene_id = columns[0]
    index = pd.read_csv(src/filename,sep=separator,usecols=[gene_id])[gene_id]
    assert index.notna().all() and index.is_unique
    genes = set(index) if genes is None else genes.intersection(index)
    assert genes
    selected = rna.loc[rna.source.eq(tag)]
    assert set(selected.libid) <= set(columns)
    matrix = pd.read_csv(src/filename,sep=separator,usecols=[gene_id,*selected.libid]).set_index(gene_id)
    assert matrix.ge(0).all().all() and np.isfinite(matrix).all().all()
    matrix = matrix.T.groupby(selected.set_index("libid").patient).sum().T
    matrices[tag] = matrix
ordered = sorted(genes)
expression = pd.concat([matrix.loc[ordered].T for matrix in matrices.values()])
assert expression.index.is_unique and set(expression.index)==set(patients)
ranks = expression.rank(axis=1,method="average",pct=True)
catalog = pd.DataFrame({"gene":ordered,"feature":[f"P_gene_{i:05d}" for i in range(len(ordered))]})
catalog.to_csv(inputs/"gene_catalog.csv",index=False)
ranks.columns = catalog.feature
features = features.join(ranks)
future["sample_id"] = future.patient+"_CP_"+future.period.astype(str)
rows = future[["sample_id","patient"]].copy()
rows["y"] = np.log1p(future.mean_cpeptide)
design = future[["sample_id","patient","period"]].join(features,on="patient",validate="many_to_one")
design["C_period2"] = design.period.eq(2).astype(int)
design["C_period3"] = design.period.eq(3).astype(int)
design["C_nlr_late"] = design.C_nlr*(design.period-1)/2
for column in ["C_baseline","C_rituximab","C_abatacept"]:
    design[column+"_t2"] = design[column]*design.C_period2
    design[column+"_t3"] = design[column]*design.C_period3
audit = future[["sample_id","patient","period","day","integral","mean_cpeptide"]].copy()
audit["baseline_day"] = audit.patient.map(baseline.day)
audit["rna_day"] = audit.patient.map(day)
audit["cbc_day"] = audit.patient.map(paired.day)
boundaries = future.integral.le(6+1e-12).groupby(future.period).mean().reindex([1,2,3])
assays.to_csv(inputs/"discovery_assays.csv",index=False)
rna[["patient","libid","source","day"]].to_csv(inputs/"used_libraries.csv",index=False)
finish_panel(study,rows,design.drop(columns=["patient","period"]),audit,
  additional_criteria={"reporting_boundary":bool(boundaries.lt(.10).all())},
  additional_summary={"reporting_boundary_by_period":boundaries.to_dict(),"reporting_floor_nmol_l":floor,
                      "common_baseline_genes":len(genes),"reserved_trials":["TN02","START"]})
