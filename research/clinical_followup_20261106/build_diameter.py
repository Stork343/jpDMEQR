from pathlib import Path

import numpy as np
import pandas as pd

from panel import finish_panel, study_config

study = study_config("ispy_diameter")
source = Path(study["source_directory"])
inputs = Path(study["input_directory"])
roles = pd.read_csv(inputs / "roles.csv", dtype={"source_id": str}).query("role == 'discovery'").set_index("source_id")
assert roles.patient.str.startswith("IS1_").all()
matrix = pd.read_csv(source / "ispy1_panel_expression.csv").set_index("gene").loc[:, roles.index].T
genes = pd.read_csv(source / "common_pam50.csv").probe.tolist()
assert set(matrix.columns) == set(genes) and len(genes) >= 40
scores = matrix.rank(axis=1, method="average", pct=True).where(matrix.notna().mean(axis=1).ge(.8), axis=0)
program = ["BIRC5", "CCNB1", "MKI67"]
proliferation = scores[program].mean(axis=1, skipna=False) * 10
visits = ["MRI LD Baseline", "MRI LD 1-3dAC", "MRI LD InterReg", "MRI LD PreSurg"]
fields = ["SUBJECTID", "age", "HR Pos", "Her2MostPos", *visits]
clinical = pd.read_excel(source / "ispy1_clinical.xlsx", sheet_name="TCIA Patient Clinical Subset",
                         usecols=fields, dtype={"SUBJECTID": str}).set_index("SUBJECTID").loc[roles.index]
clinical = clinical.apply(pd.to_numeric, errors="raise")
static = pd.DataFrame(index=roles.index)
static["C_age"] = clinical.age.where(clinical.age.between(18,110))
static["C_hr"] = clinical["HR Pos"].where(clinical["HR Pos"].isin([0,1]))
static["C_her2"] = clinical.Her2MostPos.where(clinical.Her2MostPos.isin([0,1]))
static["C_baseline"] = np.log1p(clinical[visits[0]].where(clinical[visits[0]].gt(0)))
static["C_proliferation"] = proliferation
static = static.join((scores.drop(columns=program) - 0.5).add_prefix("P_rank_"))
rows, design, audit = [], [], []
for period, field in enumerate(visits[1:], 1):
    observed = clinical.loc[clinical[field].ge(0) & static.C_baseline.notna()]
    trend = (period-1)/2
    features = static.loc[observed.index].copy()
    features["C_period2"] = int(period == 2)
    features["C_period3"] = int(period == 3)
    features["C_proliferation_late"] = features.C_proliferation * trend
    for covariate in ["C_baseline", "C_hr", "C_her2"]:
        features[covariate+"_t2"] = features[covariate] * int(period == 2)
        features[covariate+"_t3"] = features[covariate] * int(period == 3)
    other = features.filter(like="P_rank_")
    features = features.join((other * trend).rename(columns=lambda s: s.replace("P_rank_", "P_phase_")))
    features["sample_id"] = "MRI_" + roles.loc[observed.index, "patient"] + "_" + str(period)
    design.append(features)
    for source_id in observed.index:
        key = features.at[source_id, "sample_id"]
        patient = roles.at[source_id, "patient"]
        baseline = observed.at[source_id, visits[0]]
        followup = observed.at[source_id, field]
        rows.append({"sample_id": key, "patient": patient, "y": np.log1p(followup / baseline)})
        audit.append({"sample_id": key, "patient": patient, "period": period, "baseline_mm": baseline, "followup_mm": followup})
finish_panel(study, rows, pd.concat(design), audit, additional_summary={"common_pam50_genes": len(genes)})
