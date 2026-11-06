from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd

from panel import finish_panel, study_config

study = study_config("physionet_ph")
source = Path(study["source_directory"])
inputs = Path(study["input_directory"])
roles = pd.read_csv(inputs / "roles.csv")
catalog = pd.read_csv(source / "catalog.csv").query("site == 'A'")
assert set(catalog.patient) == set(roles.loc[roles.role.eq("discovery"), "patient"])
clinical_vitals = ["HR", "O2Sat", "Temp", "SBP", "MAP", "DBP", "Resp"]
gas = ["PaCO2", "HCO3", "BaseExcess"]
remaining = ["EtCO2", "FiO2", "SaO2", "AST", "BUN", "Alkalinephos", "Calcium", "Chloride",
             "Creatinine", "Bilirubin_direct", "Glucose", "Lactate", "Magnesium", "Phosphate", "Potassium",
             "Bilirubin_total", "TroponinI", "Hct", "Hgb", "PTT", "WBC", "Fibrinogen", "Platelets"]
static = ["Age", "Gender", "Unit1", "Unit2", "HospAdmTime"]
longitudinal = [*clinical_vitals, *gas, *remaining]
fields = [*static, *longitudinal, "ICULOS", "pH"]
frames = []
for batch, cells in catalog.groupby("batch"):
    with ZipFile(source / f"A_{batch:03}.zip") as archive:
        for patient in cells.patient:
            with archive.open(patient+".psv") as stream:
                frame = pd.read_csv(stream, sep="|", usecols=fields)
            assert frame.ICULOS.is_unique and frame.ICULOS.is_monotonic_increasing
            assert frame.Age.dropna().nunique() <= 1
            baseline = frame.loc[frame.ICULOS.between(1,12)]
            if not baseline.Age.ge(18).any() or not baseline.pH.between(6.5,8.0).any():
                continue
            frame = frame.loc[frame.ICULOS.between(1,48)].copy()
            frame["patient"] = patient
            frames.append(frame)
x = pd.concat(frames, ignore_index=True)
pre = x.loc[x.ICULOS.le(12)]
measured = x.loc[x.pH.between(6.5,8.0)].copy()
base = measured.loc[measured.ICULOS.le(12)]
future = measured.loc[measured.ICULOS.between(13,48)].copy()
future["period"] = (future.ICULOS-13)//12+1
future["distance"] = abs(future.ICULOS-(6+12*future.period))
future = future.sort_values(["patient", "period", "distance", "ICULOS"]).drop_duplicates(["patient", "period"])
pre = pre.loc[pre.patient.isin(future.patient)].sort_values(["patient","ICULOS"])
groups = pre.groupby("patient")
z = groups[[*static,*clinical_vitals,*gas]].last().add_prefix("C_last_")
baselines = base.sort_values(["patient","ICULOS"]).groupby("patient")
z["C_first_ph"] = baselines.pH.first()
z["C_last_ph"] = baselines.pH.last()
z["C_hours_since_ph"] = 12-baselines.ICULOS.last()
paired = pre.loc[pre.PaCO2.gt(0) & pre.Resp.gt(0)].drop_duplicates("patient", keep="last").set_index("patient")
z["C_paired_co2"] = (paired.PaCO2-40)/10
z["C_paired_resp"] = (paired.Resp-20)/5
z["C_ventilation"] = z.C_paired_co2*z.C_paired_resp
early = pre.loc[pre.ICULOS.le(6)].groupby("patient")[longitudinal].mean()
late = pre.loc[pre.ICULOS.gt(6)].groupby("patient")[longitudinal].mean()
values = groups[longitudinal]
amplitude = (values.max()-values.min()).where(values.count().ge(2))
features = pd.concat([z, groups[remaining].last().add_prefix("P_last_"), values.median().add_prefix("P_median_"),
                      (late-early).add_prefix("P_change_"), amplitude.add_prefix("P_range_")], axis=1)
features = features.loc[z.index]
future["sample_id"] = future.patient + "_pH_" + future.period.astype(str)
rows = future.rename(columns={"pH":"y"})[["sample_id", "patient", "y"]]
features = future[["sample_id","patient","period"]].join(features, on="patient", validate="many_to_one")
features["C_period2"] = features.period.eq(2).astype(int)
features["C_period3"] = features.period.eq(3).astype(int)
features["C_ph_t2"] = (features.C_last_ph-7.4)*features.C_period2
features["C_ph_t3"] = (features.C_last_ph-7.4)*features.C_period3
features = features.drop(columns=["patient","period"])
audit = future[["sample_id","patient","period","ICULOS"]].copy()
audit["baseline_hour"] = audit.patient.map(baselines.ICULOS.last())
audit["pair_hour"] = audit.patient.map(paired.ICULOS)
finish_panel(study, rows, features, audit)
