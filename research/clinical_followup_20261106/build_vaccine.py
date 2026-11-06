from pathlib import Path

import numpy as np
import pandas as pd

from panel import finish_panel, study_config

study = study_config("vaccine_cxcl10")
inputs = Path(study["input_directory"])
roles = pd.read_csv(inputs/"roles.csv", dtype={"patient":str}).query("role == 'discovery'").set_index("patient")
source = pd.read_csv(inputs/"source_assays.csv", dtype={"patient": str, "value_reported": str})
assert set(source.patient) <= set(roles.index) and source.unit_reported.eq("pg/ml").all()
values = pd.to_numeric(source.value_reported, errors="coerce")
source["value"] = values.where(np.isfinite(values) & values.ge(0))
observed = source.dropna(subset=["value"]).copy()
measure = observed.groupby(["patient","day","analyte"]).value.median().unstack("analyte")
baseline = measure.xs(0,level="day").dropna(subset=["CXL10_HUMAN"])
rank = pd.read_csv(inputs/"ranked_expression.csv.gz").set_index("sample_key")
catalog = pd.DataFrame({"gene":rank.columns,"feature":[f"P_gene_{i:05d}" for i in range(len(rank.columns))]})
pattern = ["IFIT1","IFIT3","ISG15","MX1","OAS1"]
extra = catalog.loc[~catalog.gene.isin([*pattern,"CXCL10"])].copy()
extra.to_csv(inputs/"gene_catalog.csv", index=False)
points = measure.reset_index().query("day>0").dropna(subset=["CXL10_HUMAN"])
points = points.loc[points.patient.isin(baseline.index)].copy()
points["sample_key"] = points.patient+"|"+points.day.astype(str)
points["base_key"] = points.patient+"|0"
points = points.loc[points.sample_key.isin(rank.index) & points.base_key.isin(rank.index)]
assert len(points) > 0 and points.groupby("patient").size().le(3).all()
points["sample_id"] = "FLU_" + points.patient + "_" + points.day.astype(str)
points = points.set_index("sample_id")
visit = rank.loc[points.sample_key].set_axis(points.index)
start = rank.loc[points.base_key].set_axis(points.index)
difference = visit-start
features = difference[extra.gene].set_axis(extra.feature, axis=1)
c = pd.DataFrame(index=points.index)
c["C_isg_change"] = difference[pattern].mean(axis=1,skipna=False)*10
c["C_baseline_isg"] = start[pattern].mean(axis=1,skipna=False)*10
c["C_baseline_cxcl10"] = np.log1p(points.patient.map(baseline.CXL10_HUMAN))
c["C_study269"] = points.patient.map(roles.study_accession).eq("SDY269").astype(int)
c["C_lai_vaccine"] = points.patient.map(roles.vaccine).str.startswith("LAIV").astype(int)
c["C_age"] = pd.to_numeric(points.patient.map(roles.age_reported))
c["C_female"] = points.patient.map(roles.gender).map({"Female":1,"Male":0})
c["C_log_day"] = np.log2(points.day/3)
c["C_study_day"] = c.C_study269 * c.C_log_day
c["C_baseline_day"] = c.C_baseline_cxcl10 * c.C_log_day
for marker in ["IFNG_HUMAN","IFNA2_HUMAN","CCL2_HUMAN","CCL5_HUMAN"]:
    if observed.loc[observed.analyte.eq(marker) & observed.day.eq(0),"study_accession"].nunique() == 2:
        c["C_base_"+marker] = np.log1p(points.patient.map(baseline[marker]))
rows = points[["patient"]].copy()
rows["y"] = np.log1p(points.CXL10_HUMAN)
audit = points[["patient","day","sample_key","base_key"]].copy()
audit["study"] = points.patient.map(roles.study_accession)
paired_days = audit.groupby("study").day.nunique()
finish_panel(study, rows.reset_index(), c.join(features).reset_index(), audit.reset_index(),
             additional_criteria={"repeated_days_in_each_study": bool(len(paired_days)==2 and paired_days.ge(2).all())},
             additional_summary={"paired_days_per_study":paired_days.to_dict(), "source_extra_genes":len(extra),
                                 "nonquantitative_serum_records":int(source.value.isna().sum())})
