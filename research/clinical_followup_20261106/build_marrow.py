import json
from pathlib import Path

import numpy as np
import pandas as pd

from panel import finish_panel, study_config

study = study_config("commpass_lightchain")
inputs = Path(study["input_directory"])
source = Path(study["source_directory"])
roles = pd.read_csv(inputs/"roles.csv",dtype=str).set_index("patient")
cohort = roles.loc[roles.role.eq("discovery")]
cases = json.loads((source/"discovery_clinical.json").read_text())
assert set(case["submitter_id"] for case in cases)==set(cohort.index)
kappa = "Serum Free Immunoglobulin Light Chain, Kappa"
lambda_name = "Serum Free Immunoglobulin Light Chain, Lambda"
proliferation = ["MKI67","TOP2A","BIRC5","CCNB1","UBE2C"]
clinical_units = {"Albumin":"g/L","Creatinine":"umol/L","Hemoglobin":"mmol/L","Beta 2 Microglobulin":"mcg/mL",
                  "Lactate Dehydrogenase":"ukat/L","M Protein":"g/dL"}
lab_columns = {"Albumin":"C_albumin","Creatinine":"C_creatinine","Hemoglobin":"C_hb","Beta 2 Microglobulin":"C_beta2",
               "Lactate Dehydrogenase":"C_ldh","M Protein":"C_mprotein"}
therapy_terms = ["Bortezomib","Lenalidomide","Dexamethasone","Carfilzomib","Daratumumab","Cyclophosphamide"]
tests, demographics = [], []
for case in cases:
    patient = case["submitter_id"]
    drugs = [t.get("therapeutic_agents","") for d in case.get("diagnoses",[]) for t in d.get("treatments",[])
             if t.get("days_to_treatment_start") is not None and t["days_to_treatment_start"]<=30]
    start = {"patient":patient, "C_age":case.get("demographic",{}).get("age_at_index"),
             "C_documented_drugs":int(bool(drugs))}
    for drug in therapy_terms:
        start["C_"+drug.lower()] = int(any(drug.lower() in name.lower() for name in drugs if name))
    demographics.append(start)
    for visit in case.get("follow_ups",[]):
        day = visit.get("days_to_follow_up")
        if day is None:
            continue
        for test in visit.get("molecular_tests",[]):
            name = test.get("laboratory_test")
            if name not in [kappa,lambda_name,*clinical_units] or test.get("test_value") is None:
                continue
            if test.get("test_result") != "Test Value Reported":
                continue
            amount = float(test["test_value"])
            if not np.isfinite(amount) or amount<0:
                continue
            if name in [kappa,lambda_name]:
                conversion = {"mg/dL":10,"mg/L":1,"mcg/mL":1,"g/L":1000}.get(test.get("test_units"))
                if conversion is None or test.get("biospecimen_type")!="Blood":
                    continue
                amount *= conversion
            elif test.get("test_units")!=clinical_units[name] or (name=="M Protein" and test.get("biospecimen_type")!="Blood"):
                continue
            measured = test.get("days_to_test")
            tests.append({"patient":patient,"visit":visit["follow_up_id"],"day":day,
                          "measured":day if measured is None else measured,"name":name,"value":amount})
tests = pd.DataFrame(tests)
base = tests.loc[tests.day.between(-30,30) & tests.measured.le(30)].copy()
pairs = base.loc[base.name.isin([kappa,lambda_name])].pivot_table(index=["patient","day","visit"],columns="name",values="value",aggfunc="median")
pairs = (pairs.dropna(subset=[kappa,lambda_name]).reset_index()
         .groupby(["patient","day"],as_index=False)[[kappa,lambda_name]].median().sort_values(["patient","day"]))
pairs = pairs.drop_duplicates("patient",keep="last").set_index("patient")
ratio = pairs[kappa]/pairs[lambda_name].where(pairs[lambda_name].gt(0))
pairs = pairs.loc[ratio.gt(1.65)|ratio.lt(.26)].copy()
pairs["chain"] = np.where(ratio.loc[pairs.index].gt(1.65),kappa,lambda_name)
pairs["involved"] = np.where(pairs.chain.eq(kappa),pairs[kappa],pairs[lambda_name])
pairs["uninvolved"] = np.where(pairs.chain.eq(kappa),pairs[lambda_name],pairs[kappa])
pairs = pairs.loc[pairs.involved.ge(100) & pairs.uninvolved.ge(0)]
clinical = pd.DataFrame(demographics).set_index("patient").loc[pairs.index]
clinical["C_involved"] = np.log1p(pairs.involved)
clinical["C_uninvolved"] = np.log1p(pairs.uninvolved)
clinical["C_kappa"] = pairs.chain.eq(kappa).astype(int)
last = base.groupby(["patient","name","day"],as_index=False).value.median().sort_values("day")
last = last.drop_duplicates(["patient","name"],keep="last").pivot(index="patient",columns="name",values="value")
for lab,feature in lab_columns.items():
    clinical[feature] = last[lab]
future = tests.loc[tests.patient.isin(pairs.index) & tests.name.eq(tests.patient.map(pairs.chain))
                   & tests.day.between(45,315)].copy()
periods = lambda x: np.select([x.between(45,135),x.between(136,225),x.between(226,315)],[1,2,3],default=0)
future["period"] = periods(future.day)
future = future.loc[future.period.eq(periods(future.measured)) & future.period.gt(0)]
future = future.groupby(["patient","day","period"],as_index=False).value.median()
future["distance"] = (future.day-90*future.period).abs()
future = future.sort_values(["patient","period","distance","day"]).drop_duplicates(["patient","period"])
patients = pd.Index(sorted(future.patient.unique()),name="patient")
sample = cohort.loc[patients,"sample"]
annotation = pd.read_csv(source/"first_star.tsv",sep="\t",comment="#",usecols=["gene_id","gene_name","gene_type"])
annotation = annotation.loc[annotation.gene_type.eq("protein_coding")].set_index("gene_id").gene_name
assert annotation.index.is_unique
expr = pd.read_csv(source/"transcripts.tsv.gz",sep="\t",usecols=["Ensembl_ID",*sample]).set_index("Ensembl_ID")
expr = expr.loc[expr.index.intersection(annotation.index)].groupby(annotation).median().T.loc[sample]
expr.index = patients
assert set(proliferation).issubset(expr.columns)
clinical["C_proliferation"] = expr[proliferation].mean(axis=1,skipna=False)
expr = expr.drop(columns=proliferation)
catalog = pd.DataFrame({"gene":expr.columns,"feature":[f"P_gene_{i:05d}" for i in range(len(expr.columns))]})
expr.columns = catalog.feature
catalog.to_csv(inputs/"gene_catalog.csv",index=False)
rows, audit = [], []
for v in future.itertuples(index=False):
    key = f"MM_{v.patient}_{v.period}"
    rows.append({"patient":v.patient,"sample_id":key,"y":np.log1p(v.value/pairs.at[v.patient,"involved"])})
    audit.append({"patient":v.patient,"sample_id":key,"sample":sample.at[v.patient],"period":v.period,"day":v.day,
                  "base_day":pairs.at[v.patient,"day"],"chain":pairs.at[v.patient,"chain"],
                  "baseline_mg_l":pairs.at[v.patient,"involved"],"future_mg_l":v.value})
index = pd.DataFrame(audit).set_index("sample_id")
x = index[["patient"]].join(clinical.join(expr),on="patient").drop(columns="patient")
x["C_period2"] = index.period.eq(2).astype(int)
x["C_period3"] = index.period.eq(3).astype(int)
x["C_involved_t2"] = x.C_involved*x.C_period2
x["C_involved_t3"] = x.C_involved*x.C_period3
x["C_proliferation_late"] = x.C_proliferation*(index.period-1)/2
zeros = future.value.eq(0).groupby(future.period).mean().reindex([1,2,3])
finish_panel(study,rows,x.reset_index(),audit,additional_criteria={"zero_fraction":bool(zeros.lt(.05).all())},
             additional_summary={"zero_by_period":zeros.to_dict(),"protein_coding_gene_features":len(expr.columns)})
