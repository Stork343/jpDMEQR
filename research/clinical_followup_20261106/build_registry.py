from pathlib import Path
from zipfile import ZipFile
import json

import numpy as np
import pandas as pd

from panel import study_config, finish_panel

study = study_config("t1dx_monitoring")
inputs = Path(study["input_directory"])
source = Path(study["source_directory"])
roles = pd.read_csv(inputs / "roles.csv", dtype={"patient": str})
cohort = roles.loc[roles.role.eq("discovery")].set_index("patient")
with ZipFile(source / "registry.zip") as archive:
    visit_file, = [n for n in archive.namelist() if n.endswith("/VisitsLong.txt")]
    lab_file, = [n for n in archive.namelist() if n.endswith("/Hba1cLong.txt")]
    visits = pd.read_csv(archive.open(visit_file), sep="|", dtype=str)
    labs = pd.read_csv(archive.open(lab_file), sep="|", dtype=str)
visits = visits.loc[visits.PtID.isin(cohort.index) & visits.Visit.isin(["Enrollment1","Enrollment2"])].copy()
visits["ExamDtMnC"] = pd.to_numeric(visits.ExamDtMnC)
visits = visits.merge(cohort.reset_index()[["patient","Visit","ExamDtMnC"]],
    left_on=["PtID","Visit","ExamDtMnC"], right_on=["patient","Visit","ExamDtMnC"], validate="many_to_one")
visits = visits.drop_duplicates().set_index("patient")
assert visits.index.is_unique and set(visits.index) == set(cohort.index)
landmark = visits.ExamDtMnC.clip(lower=0)
labs = labs.loc[labs.PtID.isin(cohort.index), ["PtID","HbA1c","HbA1cImputeDtMnC"]].copy()
labs["value"] = pd.to_numeric(labs.HbA1c, errors="raise")
labs["month"] = pd.to_numeric(labs.HbA1cImputeDtMnC, errors="raise")
labs = labs.dropna(subset=["value","month"]).groupby(["PtID","month"], as_index=False).value.median()
assert labs.value.gt(0).all() and labs.value.le(100).all()
labs["elapsed"] = labs.month - labs.PtID.map(landmark)
prior = labs.loc[labs.elapsed.between(-6,0)].sort_values(["PtID","month"])
last = prior.drop_duplicates("PtID", keep="last").set_index("PtID")
first = prior.drop_duplicates("PtID").set_index("PtID")
matrix = pd.DataFrame(index=visits.index)


def number(field):
    return pd.to_numeric(visits[field], errors="raise")


matrix["C_age"] = number("ageConsent")
matrix["C_duration"] = number("diabDur")
matrix["C_female"] = visits.Pt_Gender.map({"Female":1, "Male":0})
matrix["C_bmi"] = number("bmi")
matrix["C_lab"] = last.value
matrix["C_lab_change"] = last.value-first.value
matrix["C_lab_repeated"] = last.month.ne(first.month).astype(int)
matrix["C_lab_age"] = -last.elapsed
matrix["C_meter"] = number("BGTestAvgNumMeter")/2
matrix["C_check_gap"] = (number("BGTestAvgNumPtRep")-number("BGTestAvgNumMeter"))/2
matrix["C_cgm"] = visits.CGMUsed.map({"Yes":1,"No":0})
matrix["C_daily_insulin"] = number("TotalDailyInsPerKg")
matrix["C_pump"] = visits.Pt_InsulinDeliveryMethod.map({"Injections/pens":0, "Insulin pump":1,
                                                        "Pump and injections/pens":1})
income = visits.Pt_AnnualInc.where(~visits.Pt_AnnualInc.str.contains("Unknown|not to|know",case=False,na=False))
income_codes = pd.get_dummies(income, dtype=float)
income_codes.columns = [f"C_income_{i:02d}" for i in range(len(income_codes.columns))]
matrix = matrix.join(income_codes.where(income.notna(), axis=0))
numeric = ["BldPrSys","BldPrDia","WeightKg","HeightCm","bmiPerc","bmiZscore","AlbCreatRat_mggNew",
           "GFR","TSH","UnitsInsBasalOrLongAct","NumPumpBolusOrShortAct","Pt_A1cGoalLev"]
for field in numeric:
    matrix["P_"+field] = number(field)
categorical = ["Pt_PumpManuf","Pt_PumpModel","Pt_CGMManuf","Pt_CGMModel","Pt_InsPriv","Pt_InsGov",
               "Pt_InsSingleService","Pt_InsNoCoverage","Pt_MaritalStatus","Pt_raceEthNew","Pt_EduLevel",
               "Pt_ImmBioFamT1D","Pt_FatherT1D","Pt_MotherT1D","Pt_SiblingT1D","Pt_ExtBioFamT1D",
               "Pt_GrandparentT1D","AlbuminuriaStatus","ACEARB","AutonomicNeuroCl","GFRIsBelow60","Pt_SHFlg",
               "Pt_DKAFlg","Pt_SmokeCurr","Pt_A1cGoal","InjMethod","CGMUsedDevRschStudy","LipidsFastingStatus",
               *[f"InsulinUsedNew{i}" for i in range(1,7)]]
encodings = {}
extra = []
for field in categorical:
    values = visits[field].replace({"Unknown":np.nan,"Unk":np.nan,"Don't Know":np.nan})
    codes = pd.get_dummies(values, dtype=float)
    encodings[field] = list(codes.columns.astype(str))
    codes.columns = [f"P_{field}_{i:02d}" for i in range(len(codes.columns))]
    extra.append(codes.where(values.notna(),axis=0))
matrix = matrix.join(pd.concat(extra,axis=1))
for lab, factor in [("LDL",38.67),("HDL",38.67),("Triglyc",88.57),("BUN",2.801)]:
    unit = visits[lab+"Units"].str.lower().str.replace(" ","",regex=False)
    matrix["P_"+lab] = number(lab)*unit.map({"mg/dl":1,"mmol/l":factor})
for lo, hi, key in [(-24,-12,"old"),(-12,-6,"middle"),(-6,0,"recent")]:
    part = labs.loc[labs.elapsed.gt(lo) & labs.elapsed.le(hi)]
    hist = part.groupby("PtID").value.agg(["mean","min","max","std"])
    for metric in hist.columns:
        matrix[f"P_hba1c_{key}_{metric}"] = hist[metric]
future = labs.loc[labs.PtID.isin(last.index) & labs.elapsed.gt(3) & labs.elapsed.le(27)].copy()
future["period"] = np.select([future.elapsed.le(9),future.elapsed.le(15)],[1,2],default=3)
future["distance"] = (future.elapsed-future.period.map({1:6,2:12,3:24})).abs()
future = future.sort_values(["PtID","period","distance","month"]).drop_duplicates(["PtID","period"])
rows, features, audit = [], [], []
for entry in future.itertuples():
    key = f"DX_{entry.PtID}_{entry.period}"
    rows.append({"sample_id":key,"patient":entry.PtID,"y":entry.value})
    features.append(dict(matrix.loc[entry.PtID].to_dict(), sample_id=key, C_period2=int(entry.period==2),
                         C_period3=int(entry.period==3)))
    audit.append({"sample_id":key,"patient":entry.PtID,"landmark":landmark.at[entry.PtID],
                  "baseline_hba1c_month":last.at[entry.PtID,"month"],"outcome_month":entry.month,"period":entry.period})
(inputs/"categorical_levels.json").write_text(json.dumps(encodings,indent=2)+"\n")
finish_panel(study,rows,features,audit)
