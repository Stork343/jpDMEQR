from pathlib import Path

import numpy as np
import pandas as pd

from panel import finish_panel, study_config

config = study_config("vitaldb_crp")
inputs = Path(config["input_directory"])
source = Path(config["source_directory"])
roles = pd.read_csv(inputs / "roles.csv", dtype={"patient": str})
cohort = roles.loc[roles.role.eq("discovery")].set_index("caseid")
common = ["Cholecystectomy", "Distal gastrectomy", "Total gastrectomy", "Anterior resection", "Low anterior resection",
          "Right hemicolectomy", "Hemihepatectomy", "Liver segmentectomy", "Pancreaticoduodenectomy"]
published = ["caseid", "subjectid", "dis", "opname", "sex", "weight", "bmi", "asa", "emop", "preop_dm", "preop_htn",
             "preop_alb", "preop_hb", "intraop_ebl", "intraop_rbc", "intraop_crystalloid", "intraop_colloid",
             "intraop_uo", "intraop_ffp", "intraop_ppf", "intraop_mdz", "intraop_ftn", "intraop_rocu", "intraop_vecu",
             "intraop_eph", "intraop_phe", "intraop_epi", "intraop_ca", "preop_ecg", "preop_pft", "position", "ane_type"]
all_cases = pd.read_csv(source / "cases.csv.gz", usecols=published, dtype={"subjectid": str})
cases = all_cases.loc[all_cases.caseid.isin(cohort.index)].set_index("caseid")
assert cases.index.is_unique and cases.subjectid.eq(cohort.loc[cases.index, "patient"]).all()
labs = pd.read_csv(source / "labs.csv.gz", dtype=str)
labs = labs.loc[labs.caseid.isin(set(cohort.index.astype(str)))].copy()
labs["caseid"] = pd.to_numeric(labs.caseid)
labs["dt"] = pd.to_numeric(labs.dt)
labs["result"] = pd.to_numeric(labs.result)
labs = labs.dropna(subset=["dt", "result"]).groupby(["caseid", "name", "dt"], as_index=False).result.median()
labs = labs.sort_values(["caseid", "name", "dt"])
start = labs.caseid.map(cohort.opstart)
end = labs.caseid.map(cohort.landmark)
before = labs.loc[labs.dt.lt(start) & labs.dt.ge(start - 30*86400)].copy()
during = labs.loc[labs.dt.ge(start) & labs.dt.le(end)].copy()
pre = before.groupby(["caseid", "name"]).result.agg(first="first", last="last", low="min", high="max", n="size", mean="mean")
pre["range"] = (pre.high - pre.low).where(pre.n.ge(2))
pre["change"] = (pre["last"] - pre["first"]).where(pre.n.ge(2))
intra = during.groupby(["caseid", "name"]).result.agg(last="last", low="min", high="max", n="size")
intra["range"] = (intra.high - intra.low).where(intra.n.ge(2))
matrix = pd.DataFrame(index=cohort.index)
matrix["C_age"] = cohort.age
matrix["C_female"] = cases.sex.map({"F": 1, "M": 0})
matrix["C_bmi"] = cases.bmi
matrix["C_ASA"] = cases.asa
matrix["C_emergency"] = cases.emop
matrix["C_diabetes"] = cases.preop_dm
matrix["C_open"] = cohort.approach.eq("Open").astype(int)
matrix["C_robotic"] = cohort.approach.eq("Robotic").astype(int)
matrix["C_duration"] = (cohort.opend - cohort.opstart) / 3600 - 3
matrix["C_open_duration"] = matrix.C_open * matrix.C_duration
matrix["C_ebl"] = cases.intraop_ebl / 100
matrix["C_rbc"] = cases.intraop_rbc.gt(0).astype(float).where(cases.intraop_rbc.notna())
matrix["C_crystalloid"] = cases.intraop_crystalloid / cases.weight.where(cases.weight.gt(0))
matrix["C_albumin"] = cases.preop_alb
matrix["C_hb"] = cases.preop_hb
known_crp = labs.loc[labs.name.eq("crp") & labs.dt.le(end) & labs.dt.ge(start - 30*86400) & labs.result.ge(0)]
matrix["C_last_crp"] = known_crp.groupby("caseid").result.last()
for i, organ in enumerate(["Biliary/Pancreas", "Hepatic", "Stomach"]):
    matrix[f"C_organ_{i}"] = cohort.optype.eq(organ).astype(int)
for i, proc in enumerate(common):
    matrix[f"C_procedure_{i}"] = cases.opname.eq(proc).astype(int)

extra = []
for tag, group, summaries in [("prior", pre, ["first", "last", "range", "change", "mean"]),
                              ("intra", intra, ["last", "range"])]:
    panel = group[summaries].unstack()
    panel.columns = [f"P_{tag}_{name.replace('%','pct')}_{metric}" for metric, name in panel.columns]
    extra.append(panel)
for var in ["opname", "preop_ecg", "preop_pft", "position", "ane_type"]:
    codes = pd.get_dummies(cases[var], dtype=float)
    codes.columns = [f"P_{var}_{j}" for j in range(len(codes.columns))]
    extra.append(codes)
for var in [c for c in cases.columns if c.startswith("intraop_") and c not in ["intraop_ebl", "intraop_crystalloid"]]:
    extra.append(cases[[var]].rename(columns={var: "P_" + var}))
matrix = matrix.join(pd.concat(extra, axis=1))
assert matrix.columns.is_unique
future = labs.loc[labs.name.eq("crp") & labs.result.ge(0) & labs.dt.le(labs.caseid.map(cases.dis))].copy()
future["hours"] = (future.dt - future.caseid.map(cohort.landmark)) / 3600
future = future.loc[future.hours.gt(12) & future.hours.le(84)].copy()
future["period"] = np.select([future.hours.le(36), future.hours.le(60)], [1, 2], default=3)
future["distance"] = (future.hours - 24 * future.period).abs()
future = future.sort_values(["caseid", "period", "distance", "dt"]).drop_duplicates(["caseid", "period"])
rows, features, audit = [], [], []
for record in future.itertuples():
    patient = cohort.at[record.caseid, "patient"]
    key = f"VITAL_{patient}_{record.period}"
    rows.append({"sample_id": key, "patient": patient, "y": record.result})
    features.append(dict(matrix.loc[record.caseid].to_dict(), sample_id=key,
                         C_period2=int(record.period == 2), C_period3=int(record.period == 3)))
    audit.append({"sample_id": key, "patient": patient, "caseid": record.caseid, "period": record.period,
                  "landmark": cohort.at[record.caseid, "landmark"], "hours_after_landmark": record.hours,
                  "outcome_time": record.dt})
finish_panel(config, rows, features, audit)
