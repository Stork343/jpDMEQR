from pathlib import Path
import sys

import numpy as np
import pandas as pd

from panel import finish_panel, study_config
from serum import genomic_inputs, marker_series

name = sys.argv[1]
study = study_config(name)
inputs = Path(study["input_directory"])
source = Path(study["source_directory"])
folder = source/"msk_chord_2024"
cohort = pd.read_csv(inputs/"roles.csv",dtype={"patient":str}).query("role=='discovery'").set_index("patient")
assert cohort.SAMPLE_ID.is_unique
cn,mutation,events = genomic_inputs(source,cohort)
target_genes = study["marker_genes"]
if name == "chord_pancreas":
    target = (mutation.SMAD4.eq(1)|cn.SMAD4.eq(-2)).astype(float).where(mutation.SMAD4.eq(1)|cn.SMAD4.notna())
else:
    hotspot = events.HGVSp_Short.isin(["p.E380Q","p.Y537S","p.Y537N","p.Y537C","p.D538G"]) if name=="chord_breast" else (
        events.HGVSp_Short.str.match(r"^p\.[A-Z](?:12|13|61|117|146)(?:[A-Z*]|del|ins|fs)",na=False))
    positive = events.loc[events.Hugo_Symbol.isin(target_genes) & hotspot,"patient"]
    target = pd.Series(cohort.index.isin(positive).astype(float),index=cohort.index)
marker = study["marker"]
unit = "ng/ml" if marker=="cea" else "units/ml"
labs, unquantified = marker_series(folder,cohort,marker,study["report_floor"],unit)
base = labs.loc[labs.elapsed.between(-60,0)].sort_values("day").drop_duplicates("PATIENT_ID",keep="last").set_index("PATIENT_ID")
base = base.loc[base.value.ge(study["baseline_threshold"])]
history = labs.loc[labs.PATIENT_ID.isin(base.index) & labs.elapsed.between(-180,0)].copy()
prior = history.loc[history.day.le(history.PATIENT_ID.map(base.day)-14)].sort_values("day").drop_duplicates("PATIENT_ID").set_index("PATIENT_ID")
f = pd.DataFrame(index=base.index)
f["C_baseline"] = np.log1p(base.value)
f["C_days_since_test"] = -base.elapsed
f["C_prior_slope30"] = 30*(f.C_baseline-np.log1p(prior.value))/(base.day-prior.day)
f["C_genotype"] = target
f["C_specimen_years"] = (cohort.index_day-cohort.specimen_day)/365.25
f["C_metastatic_specimen"] = cohort.SAMPLE_TYPE.eq("Metastasis").astype(int)
for tag in ["IMPACT410","IMPACT468","IMPACT505"]:
    f[f"C_panel_{tag}"] = cohort.GENE_PANEL.eq(tag).astype(int)
for i,agent in enumerate(study["index_agents"]):
    f[f"C_drug_{i}"] = cohort.AGENT.str.contains(agent,regex=False).astype(int)
for i,histology in enumerate(study["histologies"][1:]):
    f[f"C_histology_{i}"] = cohort.CANCER_TYPE_DETAILED.eq(histology).astype(int)
sample = pd.read_csv(folder/"data_clinical_sample.txt",sep="\t",comment="#",dtype=str,
                     usecols=["SAMPLE_ID","TUMOR_PURITY","TMB_NONSYNONYMOUS","PRIMARY_SITE"])
sample = sample.loc[sample.SAMPLE_ID.isin(cohort.SAMPLE_ID)].set_index("SAMPLE_ID")
f["C_purity"] = pd.to_numeric(cohort.SAMPLE_ID.map(sample.TUMOR_PURITY), errors="raise")
tmb = pd.to_numeric(cohort.SAMPLE_ID.map(sample.TMB_NONSYNONYMOUS),errors="raise")
f["C_tmb"] = np.log1p(tmb.where(tmb.ge(0)))
sites = cohort.SAMPLE_ID.map(sample.PRIMARY_SITE)
for i,site in enumerate(sorted(sites.dropna().unique())):
    f[f"C_primary_site_{i}"] = sites.eq(site).astype(float).where(sites.notna())
therapy = pd.read_csv(folder/"data_timeline_treatment.txt",sep="\t",dtype=str,usecols=["PATIENT_ID","START_DATE","AGENT"])
therapy = therapy.loc[therapy.PATIENT_ID.isin(cohort.index)].copy()
earlier = therapy.loc[pd.to_numeric(therapy.START_DATE).lt(therapy.PATIENT_ID.map(cohort.index_day))]
f["C_prior_agents"] = earlier.groupby("PATIENT_ID").AGENT.nunique().reindex(f.index,fill_value=0)
other = [g for g in cn.columns if g not in target_genes]
f = f.join(cn[other].add_prefix("P_cn_")).join(mutation[other].add_prefix("P_mut_"))
future = labs.loc[labs.PATIENT_ID.isin(base.index) & labs.elapsed.between(14,112)].copy()
future["period"] = np.select([future.elapsed.le(42),future.elapsed.le(70)],[1,2],default=3)
future["distance"] = (future.elapsed-28*future.period).abs()
future = future.sort_values(["PATIENT_ID","period","distance","day"]).drop_duplicates(["PATIENT_ID","period"])
rows,features,audit = [],[],[]
for v in future.itertuples(index=False):
    pt = v.PATIENT_ID
    t = (v.period-1)/2
    key = f"{marker}_{pt}_{v.period}"
    x = dict(f.loc[pt].to_dict(),sample_id=key,C_period2=int(v.period==2),C_period3=int(v.period==3),
             C_genotype_late=target[pt]*t)
    for gene in ["TP53","RB1","PTEN","BRCA2"]:
        x["P_late_"+gene] = t*cn.at[pt,gene]
    rows.append({"sample_id":key,"patient":pt,"y":np.log1p(v.value/base.at[pt,"value"])})
    features.append(x)
    audit.append({"sample_id":key,"patient":pt,"sample":cohort.at[pt,"SAMPLE_ID"],"specimen_day":cohort.at[pt,"specimen_day"],
        "index_day":cohort.at[pt,"index_day"],"baseline_day":base.at[pt,"day"],"outcome_day":v.day,
        "period":v.period,"at_limit":v.value==study["report_floor"]})
bounds = pd.DataFrame(audit).groupby("period").at_limit.mean().reindex([1,2,3])
group_size = target.reindex(future.PATIENT_ID.unique()).value_counts().reindex([0,1],fill_value=0)
finish_panel(study,rows,features,audit,additional_criteria={"boundary":bool(bounds.lt(.05).all()),
             "genomic_groups":bool(group_size.ge(12).all())},additional_summary={"boundary_by_period":bounds.to_dict(),
             "clinical_genotype_groups":group_size.to_dict(),"common_assayed_genes":len(cn.columns),
             "marker_records_without_unique_reportable_value":unquantified})
