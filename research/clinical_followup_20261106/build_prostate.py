import json
from pathlib import Path

import numpy as np
import pandas as pd

from panel import finish_panel, study_config

study = study_config("chord_psa")
inputs = Path(study["input_directory"])
source = Path(study["source_directory"])
folder = source / "msk_chord_2024"
roles = pd.read_csv(inputs / "roles.csv", dtype={"patient":str})
cohort = roles.loc[roles.role.eq("discovery")].set_index("patient")
samples = cohort.SAMPLE_ID.to_dict()
panels = []
for name in ["IMPACT341","IMPACT410","IMPACT468","IMPACT505"]:
    record = json.loads((source/"gene_panels"/f"{name}.json").read_text())
    assert record["genePanelId"] == name
    panels.append({g["hugoGeneSymbol"] for g in record["genes"]})
genes = sorted(set.intersection(*panels))
assert {"AR","TP53","RB1","PTEN","BRCA2"}.issubset(genes)
cna_path = folder/"data_cna.txt"
cn_columns = pd.read_csv(cna_path,sep="\t",comment="#",nrows=0).columns
columns = [k for k in samples.values() if k in cn_columns]
cn = pd.read_csv(cna_path,sep="\t",comment="#",usecols=["Hugo_Symbol",*columns],index_col="Hugo_Symbol")
assert cn.index.is_unique
cn = cn.reindex(genes).T.reindex(list(samples.values()))
assert np.isin(cn.stack().to_numpy(),[-2,-1.5,-1,0,1,2]).all()
altered = pd.DataFrame(0.0,index=cn.index,columns=genes)
allowed = {"Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del",
           "In_Frame_Ins","Splice_Site","Nonstop_Mutation","Translation_Start_Site"}
for chunk in pd.read_csv(folder/"data_mutations.txt",sep="\t",comment="#",dtype=str,chunksize=50000,
                         usecols=["Hugo_Symbol","Tumor_Sample_Barcode","Variant_Classification"]):
    q = chunk.loc[chunk.Tumor_Sample_Barcode.isin(altered.index) & chunk.Hugo_Symbol.isin(genes)
                  & chunk.Variant_Classification.isin(allowed), ["Tumor_Sample_Barcode","Hugo_Symbol"]]
    for sample, gene in q.drop_duplicates().itertuples(index=False, name=None):
        altered.at[sample,gene] = 1
rename = {s:patient for patient,s in samples.items()}
cn = cn.rename(index=rename)
altered = altered.rename(index=rename)
tp53 = (altered.TP53.eq(1) | cn.TP53.eq(-2)).astype(float).where(altered.TP53.eq(1) | cn.TP53.notna())
ar = cn.AR.ge(2).astype(float).where(cn.AR.notna())
psa = pd.read_csv(folder/"data_timeline_psa_labs.txt",sep="\t",comment="#",dtype=str,
                  usecols=["PATIENT_ID","START_DATE","TEST","RESULT","LR_UNIT_MEASURE"])
psa = psa.loc[psa.PATIENT_ID.isin(cohort.index) & psa.TEST.eq("PSA")].copy()
assert psa.LR_UNIT_MEASURE.str.strip().str.lower().eq("ng/ml").all()
psa["day"] = pd.to_numeric(psa.START_DATE,errors="raise")
groups = psa.RESULT.str.extract(r"^\s*(<=|<|≤|=|>=|>|≥)?\s*(\d+(?:\.\d*)?|\.\d+)\s*$")
values = pd.to_numeric(groups[1])
exact = groups[0].isna() | groups[0].eq("=")
low = groups[0].isin(["<","<=","≤"]) & values.le(.1)
psa["reportable"] = values.clip(lower=.1).where(exact | low)
psa = psa.dropna(subset=["day","reportable"]).groupby(["PATIENT_ID","day"],as_index=False).reportable.median()
psa["elapsed"] = psa.day-psa.PATIENT_ID.map(cohort.index_day)
baseline = psa.loc[psa.elapsed.between(-60,0)].sort_values("day").drop_duplicates("PATIENT_ID",keep="last").set_index("PATIENT_ID")
baseline = baseline.loc[baseline.reportable.ge(2)]
earlier = psa.loc[psa.elapsed.between(-180,0) & psa.PATIENT_ID.isin(baseline.index)].copy()
earlier = earlier.loc[earlier.day.le(earlier.PATIENT_ID.map(baseline.day)-14)].sort_values("day").drop_duplicates("PATIENT_ID").set_index("PATIENT_ID")
features = pd.DataFrame(index=baseline.index)
features["C_psa"] = np.log1p(baseline.reportable)
features["C_psa_age"] = -baseline.elapsed
features["C_psa_slope30"] = (features.C_psa-np.log1p(earlier.reportable))*30/(baseline.day-earlier.day)
features["C_tp53"] = tp53
features["C_ar_amp"] = ar
features["C_met_sample"] = cohort.SAMPLE_TYPE.eq("Metastasis").astype(int)
features["C_specimen_age"] = (cohort.index_day-cohort.specimen_day)/365.25
for agent in ["ENZALUTAMIDE","APALUTAMIDE","DAROLUTAMIDE"]:
    features["C_"+agent.lower()] = cohort.AGENT.eq(agent).astype(int)
for name in ["IMPACT410","IMPACT468","IMPACT505"]:
    features["C_"+name] = cohort.GENE_PANEL.eq(name).astype(int)
treatment = pd.read_csv(folder/"data_timeline_treatment.txt",sep="\t",dtype=str,usecols=["PATIENT_ID","START_DATE","AGENT"])
treatment = treatment.loc[treatment.PATIENT_ID.isin(cohort.index)].copy()
prior_taxane = treatment.loc[pd.to_numeric(treatment.START_DATE).lt(treatment.PATIENT_ID.map(cohort.index_day))
                            & treatment.AGENT.str.contains(r"DOCETAXEL|CABAZITAXEL",case=False,na=False)]
features["C_prior_taxane"] = features.index.isin(prior_taxane.PATIENT_ID).astype(int)
grade = pd.read_csv(folder/"data_timeline_gleason.txt",sep="\t",dtype=str,usecols=["PATIENT_ID","START_DATE","GLEASON_SCORE"])
grade = grade.loc[grade.PATIENT_ID.isin(cohort.index)].copy()
grade["day"] = pd.to_numeric(grade.START_DATE)
grade["score"] = pd.to_numeric(grade.GLEASON_SCORE,errors="raise")
grade = grade.loc[grade.day.le(grade.PATIENT_ID.map(cohort.index_day)) & grade.score.between(2,10)]
grade = grade.groupby(["PATIENT_ID","day"],as_index=False).score.median().sort_values("day").drop_duplicates("PATIENT_ID",keep="last")
features["C_gleason"] = grade.set_index("PATIENT_ID").score
other = [g for g in genes if g!="TP53"]
features = features.join(cn[other].add_prefix("P_cn_")).join(altered[other].add_prefix("P_mut_"))
future = psa.loc[psa.PATIENT_ID.isin(baseline.index) & psa.elapsed.between(14,112)].copy()
future["period"] = np.select([future.elapsed.le(42),future.elapsed.le(70)],[1,2],default=3)
future["distance"] = (future.elapsed-28*future.period).abs()
future = future.sort_values(["PATIENT_ID","period","distance","day"]).drop_duplicates(["PATIENT_ID","period"])
rows, matrix, audit = [], [], []
for visit in future.itertuples():
    patient = visit.PATIENT_ID
    trend = (visit.period-1)/2
    key = f"PSA_{patient}_{visit.period}"
    x = dict(features.loc[patient].to_dict(),sample_id=key,C_period2=int(visit.period==2),C_period3=int(visit.period==3),
             C_tp53_late=features.at[patient,"C_tp53"]*trend)
    for gene in ["AR","RB1","PTEN","BRCA2"]:
        x["P_time_cn_"+gene] = cn.at[patient,gene]*trend
    matrix.append(x)
    rows.append({"sample_id":key,"patient":patient,"y":np.log1p(visit.reportable/baseline.at[patient,"reportable"])})
    audit.append({"sample_id":key,"patient":patient,"sample":cohort.at[patient,"SAMPLE_ID"],"period":visit.period,
                  "sample_day":cohort.at[patient,"specimen_day"],"index_day":cohort.at[patient,"index_day"],
                  "baseline_day":baseline.at[patient,"day"],"outcome_day":visit.day,"at_limit":visit.reportable==.1})
boundaries = pd.DataFrame(audit).groupby("period").at_limit.mean().reindex([1,2,3])
groups = tp53.reindex(future.PATIENT_ID.unique()).value_counts().reindex([0,1],fill_value=0)
finish_panel(study,rows,matrix,audit,
             additional_criteria={"lower_boundary":bool(boundaries.lt(.05).all()),"genomic_groups":bool(groups.ge(12).all())},
             additional_summary={"lower_boundary_fractions":boundaries.to_dict(),"tp53_groups":groups.to_dict(),
                                 "common_assayed_genes":len(genes)})
