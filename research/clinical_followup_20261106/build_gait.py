from pathlib import Path

import numpy as np
import pandas as pd

from panel import finish_panel, study_config

config = study_config("gaitrec_propulsion")
inputs = Path(config["input_directory"])
source = Path(config["source_directory"])
roles = pd.read_csv(inputs / "roles.csv", dtype={"patient": str})
ids = set(roles.loc[roles.role.eq("discovery"), "patient"])
meta = pd.read_csv(source / "GRF_metadata.csv", dtype={"SUBJECT_ID": str})
meta = meta.loc[meta.SUBJECT_ID.isin(ids) & meta.TRAIN.eq(1) & meta.TEST.eq(0) & meta.SPEED.eq(2)
                & meta.READMISSION.eq(0) & meta.AGE.ge(18) & meta.AFFECTED_SIDE.isin([0, 1])
                & meta.SHOD_CONDITION.notna() & meta.ORTHOPEDIC_INSOLE.notna()].copy()
meta["date"] = pd.to_datetime(meta.SESSION_DATE, format="%d.%m.%Y")
first_date = meta.loc[meta.SESSION_TYPE.eq(1)].groupby("SUBJECT_ID").date.min()
meta["day"] = (meta.date - meta.SUBJECT_ID.map(first_date)).dt.days
base_sessions = meta.loc[meta.SESSION_TYPE.eq(1) & meta.day.eq(0)]
follow = meta.loc[meta.SESSION_TYPE.eq(2) & meta.day.between(4, 35)].copy()
assert meta.SESSION_ID.is_unique
keys = ["SUBJECT_ID", "SESSION_ID", "TRIAL_ID"]
signals = ["F_AP", "F_V", "F_ML", "COP_AP", "COP_ML"]
curves = {}


def load_curve(signal, side):
    path = source / f"GRF_{signal}_PRO_{side}.csv"
    assert path.with_name(path.name + ".verified.json").exists()
    allowed = set(base_sessions.SESSION_ID.astype(str))
    if signal == "F_AP":
        allowed.update(follow.SESSION_ID.astype(str))
    records = []
    for chunk in pd.read_csv(path, dtype=str, chunksize=5000):
        frame = chunk.loc[chunk.SUBJECT_ID.isin(ids) & chunk.SESSION_ID.isin(allowed)].copy()
        if frame.empty:
            continue
        frame["SESSION_ID"] = pd.to_numeric(frame.SESSION_ID)
        frame["TRIAL_ID"] = pd.to_numeric(frame.TRIAL_ID)
        channels = [f"{signal}_PRO_{i}" for i in range(1, 102)]
        assert set(frame.columns) == set(keys + channels)
        values = frame.set_index(keys)[channels].apply(pd.to_numeric, errors="raise").dropna()
        values.columns = np.arange(101)
        records.append(values)
    curve = pd.concat(records).sort_index()
    assert curve.index.is_unique and np.isfinite(curve.to_numpy()).all()
    return curve


for signal in signals:
    for side in ["left", "right"]:
        curves[signal, side] = load_curve(signal, side)
baseline_trials = curves["F_V", "left"].index
for curve in curves.values():
    baseline_trials = baseline_trials.intersection(curve.index)
baseline_trials = baseline_trials[baseline_trials.get_level_values("SESSION_ID").isin(base_sessions.SESSION_ID)]
ap_trials = curves["F_AP", "left"].index.intersection(curves["F_AP", "right"].index)
base_counts = pd.Series(1, index=baseline_trials).groupby(level="SESSION_ID").size()
follow_counts = pd.Series(1, index=ap_trials).groupby(level="SESSION_ID").size()
base = base_sessions.loc[base_sessions.SESSION_ID.map(base_counts).ge(4)].sort_values("SESSION_ID").drop_duplicates("SUBJECT_ID")
follow = follow.loc[follow.SESSION_ID.map(follow_counts).ge(4)].copy()
base = base.set_index("SUBJECT_ID")
for col in ["AFFECTED_SIDE", "SHOD_CONDITION", "ORTHOPEDIC_INSOLE"]:
    follow = follow.loc[follow[col].eq(follow.SUBJECT_ID.map(base[col]))]
follow["period"] = np.select([follow.day.le(11), follow.day.le(21)], [1, 2], default=3)
follow["distance"] = (follow.day - follow.period.map({1: 7, 2: 14, 3: 28})).abs()
follow = follow.sort_values(["SUBJECT_ID", "period", "distance", "day", "SESSION_ID"]).drop_duplicates(["SUBJECT_ID", "period"])
baseline_profiles = {key: curve.loc[baseline_trials] for key, curve in curves.items()}
detailed = sorted(base.CLASS_LABEL_DETAILED.dropna().unique())
phase = np.arange(30, 96)
rows, features, audit = [], [], []
for patient, future in follow.groupby("SUBJECT_ID", sort=True):
    b = base.loc[patient]
    affected = "left" if b.AFFECTED_SIDE == 0 else "right"
    other = "right" if affected == "left" else "left"
    f = {"C_age": b.AGE, "C_male": b.SEX, "C_bmi": b.BODY_MASS / (b.HEIGHT / 100) ** 2 if b.HEIGHT > 0 else np.nan,
         "C_orthopedic_shoe": int(b.SHOD_CONDITION == 2), "C_barefoot": int(b.SHOD_CONDITION == 0),
         "C_insole": b.ORTHOPEDIC_INSOLE}
    for level in sorted(set(base.CLASS_LABEL) - {"A"}):
        f[f"C_site_{level}"] = int(b.CLASS_LABEL == level)
    for level in detailed:
        f[f"C_diagnosis_{level}"] = int(b.CLASS_LABEL_DETAILED == level) if pd.notna(b.CLASS_LABEL_DETAILED) else np.nan
    centroids = {}
    for tag, side in [("affected", affected), ("other", other)]:
        profile = {s: baseline_profiles[s, side].loc[(patient, b.SESSION_ID)] for s in signals}
        ap, vertical, lateral = profile["F_AP"], profile["F_V"], profile["F_ML"]
        positive = ap.loc[:, phase].clip(lower=0)
        mass = positive.sum(axis=1)
        centroids[tag] = positive.mul(phase, axis=1).sum(axis=1) / mass.where(mass.gt(0))
        f[f"C_{tag}_propulsion"] = 100 * ap.loc[:, 50:95].max(axis=1).clip(lower=0).median()
        f[f"C_{tag}_braking"] = 100 * (-ap.loc[:, 5:50].min(axis=1)).clip(lower=0).median()
        f[f"C_{tag}_mean_forward"] = 100 * ap.clip(lower=0).mean(axis=1).median()
        f[f"C_{tag}_vertical_peak"] = 100 * vertical.loc[:, 5:95].max(axis=1).median()
        f[f"C_{tag}_vertical_mid"] = 100 * vertical.loc[:, 35:65].min(axis=1).median()
        f[f"C_{tag}_lateral"] = 100 * lateral.loc[:, 5:95].abs().max(axis=1).median()
        for s, curve in profile.items():
            points = curve.loc[:, np.arange(4, 97, 4)].median(axis=0)
            f.update({f"P_{tag}_{s}_{int(t):02d}": value for t, value in points.items()})
    f["C_propulsion_delay"] = (centroids["affected"] - centroids["other"]).median() / 10
    for visit in future.itertuples():
        trial_ids = ap_trials[ap_trials.get_level_values("SESSION_ID") == visit.SESSION_ID]
        peak = curves["F_AP", affected].loc[trial_ids, 50:95].max(axis=1).clip(lower=0)
        key = f"GAIT_{patient}_{visit.period}"
        rows.append({"sample_id": key, "patient": patient, "y": float(100 * peak.median())})
        features.append(dict(f, sample_id=key, C_period2=int(visit.period == 2), C_period3=int(visit.period == 3)))
        audit.append({"sample_id": key, "patient": patient, "baseline_session": b.SESSION_ID,
                      "follow_session": visit.SESSION_ID, "affected": affected, "day": int(visit.day),
                      "period": int(visit.period), "trials": len(peak)})
finish_panel(config, rows, features, audit)
