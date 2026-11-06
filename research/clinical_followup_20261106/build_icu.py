from pathlib import Path
import re

import numpy as np
import pandas as pd

from panel import finish_panel, study_config

studies = {k: study_config(k) for k in ["eicu_lactate", "eicu_creatinine"]}
base = studies["eicu_lactate"]
source = Path(base["source_directory"])
roles = pd.read_csv(Path(base["input_directory"]) / "roles.csv", dtype={"patient": str})
assert roles.equals(pd.read_csv(Path(studies["eicu_creatinine"]["input_directory"]) / "roles.csv", dtype={"patient": str}))
roster = roles.loc[roles.role.eq("discovery")].set_index("patientunitstayid")
stay = "patientunitstayid"
stay_end = pd.read_csv(source / "patient.csv.gz", usecols=[stay, "unitdischargeoffset"]).set_index(stay).unitdischargeoffset
units = {"albumin": "g/dL", "ALT (SGPT)": "Units/L", "AST (SGOT)": "Units/L", "HCO3": "mmol/L",
         "bicarbonate": "mmol/L", "anion gap": "", "BUN": "mg/dL", "calcium": "mg/dL", "chloride": "mmol/L",
         "creatinine": "mg/dL", "glucose": "mg/dL", "Hct": "%", "Hgb": "g/dL", "potassium": "mmol/L",
         "sodium": "mmol/L", "WBC x 1000": "K/mcL", "platelets x 1000": "K/mcL", "MCV": "fL", "RDW": "%",
         "RBC": "M/mcL", "magnesium": "mg/dL", "phosphate": "mg/dL", "total bilirubin": "mg/dL",
         "PT - INR": "ratio", "PTT": "sec", "total protein": "g/dL", "pH": "", "paO2": "mm Hg",
         "paCO2": "mm Hg", "O2 Sat (%)": "%", "Base Excess": "mEq/L", "lactate": "mmol/L"}
vital_ranges = {"temperature": (20, 45), "sao2": (0, 100), "heartrate": (0, 300), "respiration": (0, 100),
                "cvp": (-20, 100), "etco2": (0, 180), "systemicsystolic": (0, 350),
                "systemicdiastolic": (0, 300), "systemicmean": (0, 350), "noninvasivesystolic": (0, 350),
                "noninvasivediastolic": (0, 300), "noninvasivemean": (0, 350)}


def load_table(name, columns):
    frame = pd.read_csv(source / f"{name}.csv.gz", usecols=[stay, *columns], dtype=str)
    selected = frame[stay].isin(set(roster.index.astype(str)))
    frame = frame.loc[selected].copy()
    frame[stay] = frame[stay].astype(int)
    return frame


labs = load_table("lab", ["labname", "labresultoffset", "labresultrevisedoffset", "labresult", "labresulttext",
                          "labmeasurenamesystem"])
labs = labs.loc[labs.labname.isin(units)].copy()
labs = labs.loc[labs.labmeasurenamesystem.fillna("").eq(labs.labname.map(units))
                & ~labs.labresulttext.str.contains(r"[<>≤≥]", na=False)].copy()
labs["time"] = pd.to_numeric(labs.labresultoffset)
labs["available"] = pd.to_numeric(labs.labresultrevisedoffset)
labs["value"] = pd.to_numeric(labs.labresult, errors="raise")
labs = labs.loc[labs.value.notna() & labs.time.between(0, 2880) & labs.time.le(labs[stay].map(stay_end))]
labs = labs.loc[~labs.labname.isin(["creatinine", "lactate"]) | labs.value.gt(0)]
history = labs.loc[labs.time.le(720) & labs.available.le(720), [stay, "time", "labname", "value"]].copy()
history = history.rename(columns={"labname": "variable"})
pressure = []
series = [history]
for table, variables in [("vitalPeriodic", list(vital_ranges)[:9]), ("vitalAperiodic", list(vital_ranges)[9:])]:
    frame = load_table(table, ["observationoffset", *variables])
    frame["time"] = pd.to_numeric(frame.observationoffset)
    frame = frame.loc[frame.time.between(0, 720) & frame.time.le(frame[stay].map(stay_end))].copy()
    for var in variables:
        values = pd.to_numeric(frame[var])
        low, high = vital_ranges[var]
        frame[var] = values.where(values.gt(low) & values.le(high))
    frame["hour"] = ((frame.time - 1) // 60).clip(lower=0)
    pressure.append(frame.groupby([stay, "hour"])[variables].median())
    long = frame.melt(id_vars=[stay, "time"], value_vars=variables, var_name="variable", value_name="value").dropna()
    series.append(long)

hourly = pd.concat(pressure, axis=1)
for group in ["systemic", "noninvasive"]:
    systolic = hourly[group + "systolic"]
    diastolic = hourly[group + "diastolic"]
    eligible = systolic.ge(diastolic) & diastolic.gt(0)
    hourly.loc[~eligible, [group + k for k in ["systolic", "diastolic", "mean"]]] = np.nan
chosen = {}
for field in ["diastolic", "systolic", "mean"]:
    chosen[field] = hourly["systemic" + field].where(hourly.systemicdiastolic.notna(), hourly["noninvasive" + field])
hourly["dsi"] = hourly.heartrate / chosen["diastolic"]
hourly["pulse"] = (chosen["systolic"] - chosen["diastolic"]) / 10
hourly["map"] = chosen["mean"]
signals = hourly.groupby(level=stay)[["dsi", "pulse", "map", "heartrate"]]
hemodynamics = signals.median()
hemodynamics[["dsi", "pulse"]] = hemodynamics[["dsi", "pulse"]].where(signals.count()[["dsi", "pulse"]].ge(2))

long = pd.concat(series, ignore_index=True).groupby([stay, "variable", "time"], as_index=False).value.median()
long = long.sort_values([stay, "variable", "time"])
group = long.groupby([stay, "variable"])
summaries = group.value.agg(first="first", last="last", lo="min", hi="max", n="size")
summaries["change"] = (summaries["last"] - summaries["first"]).where(summaries.n.ge(2))
summaries["range"] = (summaries.hi - summaries.lo).where(summaries.n.ge(2))
for label, mask in [("early", long.time.le(360)), ("late", long.time.gt(360))]:
    summaries[label] = long.loc[mask].groupby([stay, "variable"]).value.mean()
last_values = summaries["last"].unstack()
wide = summaries[["first", "last", "range", "change", "early", "late"]].unstack()
wide.columns = [f"P_{re.sub('[^a-z0-9]+', '_', var.lower()).strip('_')}_{stat}" for stat, var in wide.columns]
assert wide.columns.is_unique
for name, study in studies.items():
    target = "lactate" if name == "eicu_lactate" else "creatinine"
    target_summary = summaries.xs(target, level="variable")
    eligible = target_summary.index[target_summary.hi.ge(2) if target == "lactate" else target_summary.hi.gt(0)]
    future = labs.loc[labs[stay].isin(eligible) & labs.labname.eq(target) & labs.time.gt(720)].copy()
    future = future.groupby([stay, "time"], as_index=False).value.median()
    future["period"] = (future.time - 721) // 720 + 1
    future["distance"] = (future.time - (future.period * 720 + 360)).abs()
    future = future.sort_values([stay, "period", "distance", "time"]).drop_duplicates([stay, "period"])
    baseline = roster.loc[eligible, ["age", "gender"]].copy()
    baseline.columns = ["C_age", "gender"]
    baseline["C_female"] = baseline.pop("gender").map({"Female": 1, "Male": 0})
    baseline["C_first"] = target_summary["first"]
    baseline["C_last"] = target_summary["last"]
    baseline["C_change"] = target_summary["last"] - target_summary["first"]
    baseline["C_repeat"] = target_summary.n.gt(1).astype(int)
    baseline["C_map"] = hemodynamics["map"]
    baseline["C_hr"] = hemodynamics.heartrate
    if target == "lactate":
        baseline["C_dsi"] = hemodynamics.dsi
        baseline["C_peak"] = target_summary.hi
        clinical = {"C_pH": "pH", "C_creatinine": "creatinine"}
    else:
        baseline["C_pulse_pressure"] = hemodynamics.pulse
        clinical = {"C_bun": "BUN", "C_potassium": "potassium", "C_bicarbonate": "bicarbonate", "C_hgb": "Hgb"}
    for col, lab in clinical.items():
        baseline[col] = last_values[lab]
    prefix = "P_" + target + "_"
    baseline = baseline.join(wide.loc[:, ~wide.columns.str.startswith(prefix)])
    rows, features, audit = [], [], []
    for record in future.itertuples(index=False):
        stay_id = record.patientunitstayid
        patient = str(roster.at[stay_id, "patient"])
        period = int(record.period)
        key = f"{name}_{patient}_{period}"
        rows.append({"sample_id": key, "patient": patient, "y": record.value})
        features.append(dict(baseline.loc[stay_id].to_dict(), sample_id=key, C_period2=int(period == 2), C_period3=int(period == 3)))
        audit.append({"sample_id": key, "patient": patient, "stay": stay_id, "period": period, "outcome_minute": record.time})
    finish_panel(study, rows, features, audit)
