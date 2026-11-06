from pathlib import Path
import json
import sys
from zipfile import ZipFile

import pandas as pd
import yaml

config = yaml.safe_load(Path(sys.argv[1] if len(sys.argv) > 1 else "config/application_iobp2_20261106.yml").read_text())
inputs = Path(config["input_directory"])
report = Path(config["report_directory"])
roles = pd.read_csv(inputs / "roles.csv", dtype={"PtID": str})
discovery = roles.loc[roles.role.eq("discovery")].set_index("PtID")
assert len(discovery) == json.loads((report / "roles.json").read_text())["discovery_patients"]
stamp = "%m/%d/%Y %I:%M:%S %p"
randomization = pd.to_datetime(discovery.RandDt.str.split().str[0], format="%m/%d/%Y")
frames = []
counts = {"discovery_CGM_mg_dl_records": 0, "window_records": 0, "lower_coded_records": 0,
          "upper_coded_records": 0, "missing_window_glucose": 0, "missing_device_clock": 0}
columns = ["PtID", "DeviceDtTm", "RecordType", "Value", "Units"]
with ZipFile(Path(config["source_directory"]) / "IOBP2_RCT.zip") as archive:
    with archive.open("Data Tables/IOBP2DeviceCGM.txt") as stream:
        for chunk in pd.read_csv(stream, sep="|", usecols=columns, dtype=str, chunksize=300_000):
            chunk = chunk.loc[chunk.PtID.isin(discovery.index) & chunk.RecordType.eq("CGM")].copy()
            if chunk.empty:
                continue
            assert chunk.Units.eq("mg/dL").all()
            counts["discovery_CGM_mg_dl_records"] += len(chunk)
            clock = chunk.DeviceDtTm.str.contains(":", regex=False, na=False)
            counts["missing_device_clock"] += int((~clock).sum())
            chunk = chunk.loc[clock]
            times = pd.to_datetime(chunk.DeviceDtTm, format=stamp)
            assert times.notna().all()
            day = (times.dt.normalize() - chunk.PtID.map(randomization)).dt.days
            slot = (times.dt.hour * 60 + times.dt.minute) // 5
            if "forecast_day" in config:
                landmark = config["forecast_day"]
                day = day - landmark - day.le(landmark).astype(int)
            baseline = day.lt(0) if config.get("baseline_window") == "run_in" else day.between(-14, -1)
            keep = baseline | (day.between(1, config["period_days"] * config["periods"]) & slot.lt(72))
            if not keep.any():
                continue
            sample = chunk.loc[keep, ["PtID"]].copy()
            values = pd.to_numeric(chunk.loc[keep, "Value"], errors="raise")
            assert values.dropna().between(39, 401).all(), values.describe()
            counts["window_records"] += len(sample)
            counts["missing_window_glucose"] += int(values.isna().sum())
            counts["lower_coded_records"] += int(values.eq(39).sum())
            counts["upper_coded_records"] += int(values.eq(401).sum())
            sample["day"] = day.loc[keep].astype("int16")
            sample["slot"] = slot.loc[keep].astype("int16")
            sample["glucose"] = values.clip(40, 400)
            frames.append(sample.dropna())

records = pd.concat(frames, ignore_index=True)
group = records.groupby(["PtID", "day", "slot"], sort=True).glucose
slots = group.median().reset_index()
assert set(slots.PtID).issubset(discovery.index) and slots.day.ne(0).all()
assert slots.loc[slots.day.gt(0), "slot"].lt(72).all()
counts.update({"valid_window_records": len(records), "unique_five_minute_slots": len(slots),
               "extra_records_in_occupied_slots": len(records) - len(slots),
               "discovery_patients_with_window_readings": slots.PtID.nunique(),
               "reserved_patients_processed": 0, "source_table": "IOBP2DeviceCGM"})
slots.to_parquet(inputs / "discovery_slots.parquet", index=False)
(report / "source_counts.json").write_text(json.dumps(counts, indent=2) + "\n")
print(json.dumps(counts))
