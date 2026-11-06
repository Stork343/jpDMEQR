from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd

from panel import finish_panel, study_config

config = study_config("uwhvf_central")
inputs = Path(config["input_directory"])
roles = pd.read_csv(inputs / "roles.csv", dtype={"patient": str})
ids = set(roles.loc[roles.role.eq("discovery"), "patient"])
with ZipFile(Path(config["source_directory"]) / "release.zip") as archive:
    member, = [n for n in archive.namelist() if n.endswith("/CSV/Coord_242.csv")]
    coords = pd.read_csv(archive.open(member)).set_index("LocID").sort_index()
    member, = [n for n in archive.namelist() if n.endswith("/CSV/VF_Data.csv")]
    cols = ["PatID", "Gender", "Eye", "Age"] + [f"{p}_{i}" for p in ["Sens", "TD"] for i in coords.index]
    source = pd.read_csv(archive.open(member), usecols=cols, dtype=str)
    source = source.loc[source.PatID.isin(ids)].copy()
assert len(coords) == 54
blind = (coords.X.eq(15) & coords.Y.abs().eq(3))
center = coords.X.abs().eq(3) & coords.Y.abs().eq(3)
ring = ((coords.X.abs().eq(3) & coords.Y.abs().eq(9)) | (coords.X.abs().eq(9) & coords.Y.abs().eq(3)))
valid = coords.index[~blind].tolist()
peripheral = coords.index[~(blind | center)].tolist()
central_cols = [f"Sens_{i}" for i in coords.index[center]]
ring_cols = [f"TD_{i}" for i in coords.index[ring]]
td_cols = [f"TD_{i}" for i in valid]
assert len(central_cols) == 4 and len(ring_cols) == 8 and len(peripheral) == 48
numbers = [*central_cols, *td_cols]
source[numbers] = source[numbers].apply(pd.to_numeric, errors="raise").replace(100, np.nan)
source["Age"] = pd.to_numeric(source.Age, errors="raise")
assert source[central_cols].stack().between(0, 60).all()
sex = source.dropna(subset=["Gender"]).groupby("PatID").Gender.first()
exams = source.groupby(["PatID", "Eye", "Age"])[numbers].median().reset_index()
exams["center"] = exams[central_cols].mean(axis=1, skipna=False)
exams["mtd"] = exams[td_cols].mean(axis=1)
exams["baseline_ok"] = exams[td_cols].notna().mean(axis=1).ge(0.8) & exams.center.notna()
rows, features, audit = [], [], []
for patient, history in exams.groupby("PatID", sort=True):
    for eye, series in history.groupby("Eye", sort=True):
        adult = series.loc[series.Age.ge(18)].sort_values("Age")
        if adult.empty:
            continue
        first = adult.iloc[0]
        possible = adult.loc[(adult.Age - first.Age).between(0.25, 2)]
        if possible.empty:
            continue
        last = possible.iloc[0]
        if not (first.baseline_ok and last.baseline_ok and last.mtd <= -2 and last.center >= 15):
            continue
        delay = adult.Age - last.Age
        late = adult.loc[adult.center.notna() & delay.between(0.5, 4)].copy()
        if late.empty:
            continue
        fellow = history.loc[history.Eye.ne(eye) & history.Age.between(last.Age - 2, last.Age)].sort_values("Age")
        f = {"C_age": last.Age, "C_right": int(eye == "Right"),
             "C_female": {"M": 0, "F": 1}.get(sex.get(patient), np.nan), "C_interval": last.Age - first.Age,
             "C_center": last.center, "C_center_change": last.center - first.center,
             "C_field": last.mtd, "C_field_change": last.mtd - first.mtd,
             "C_field_sd": last[td_cols].astype(float).std(), "C_ring": last[ring_cols].astype(float).mean(),
             "C_fellow_center": fellow.iloc[-1].center if len(fellow) else np.nan,
             "C_fellow_field": fellow.iloc[-1].mtd if len(fellow) else np.nan}
        paired_ring = pd.concat([first[ring_cols], last[ring_cols]], axis=1).astype(float)
        f["C_ring_persistent"] = paired_ring.le(-5).all(axis=1).sum() / 2 if paired_ring.notna().all().all() else np.nan
        for tag, visit in [("index_first", first), ("index_last", last),
                           ("fellow_last", fellow.iloc[-1] if len(fellow) else None),
                           ("fellow_previous", fellow.iloc[-2] if len(fellow) >= 2 else None)]:
            for i in peripheral:
                f[f"P_{tag}_{i}"] = visit[f"TD_{i}"] if visit is not None else np.nan
        late["delay"] = late.Age - last.Age
        late["period"] = np.select([late.delay.lt(1.5), late.delay.lt(2.5)], [1, 2], default=3)
        for period, candidates in late.groupby("period"):
            selected = candidates.assign(distance=(candidates.delay - period).abs()).sort_values(["distance", "Age"]).iloc[0]
            sample = f"UW_{patient}_{eye}_{period}"
            rows.append({"sample_id": sample, "patient": patient, "y": selected.center})
            features.append(dict(f, sample_id=sample, C_period2=int(period == 2), C_period3=int(period == 3)))
            audit.append({"sample_id": sample, "patient": patient, "eye": eye, "first_age": first.Age,
                          "landmark": last.Age, "outcome_age": selected.Age, "period": int(period)})
finish_panel(config, rows, features, audit)
