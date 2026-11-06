import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd
import yaml

config = yaml.safe_load(Path("config/application_loop_20261106.yml").read_text())
source = Path(config["source_directory"])
inputs = Path(config["input_directory"])
report = Path(config["report_directory"])
assert json.loads((source / "verified.json").read_text())["all_members_verified"]
roles = pd.read_csv(inputs / "roles.csv", dtype={"PtID": str})
discovery = roles.loc[roles.role.eq("discovery")].set_index("PtID")
assert len(discovery) == json.loads((report / "roles.json").read_text())["discovery_patients"]
dates = pd.to_datetime(discovery.EnrollDt, format="%m/%d/%Y")
offsets = pd.to_timedelta(discovery.PtTimezoneOffset, unit="h")
landmark = config["forecast_day"]
end = landmark + config["periods"] * config["period_days"]
parts = inputs / "device_parts"
parts.mkdir(exist_ok=True)
columns = ["PtID", "UTCDtTm", "RecordType", "CGMVal", "Units"]
with ZipFile(source / "Loop.zip") as archive:
    assert (source / "PtRoster.txt").read_bytes() == archive.read("Data Tables/PtRoster.txt")
    for index in range(1, 7):
        cache = parts / f"device_{index}.parquet"
        metadata = cache.with_suffix(".json")
        if metadata.exists():
            assert cache.exists()
            continue
        counts = {"selected_patient_cgm": 0, "without_standard_hour": 0, "window_readings": 0,
                  "invalid_or_missing": 0, "below_range": 0, "above_range": 0}
        frames = []
        with archive.open(f"Data Tables/LOOPDeviceCGM{index}.txt") as stream:
            for chunk in pd.read_csv(stream, sep="|", dtype=str, usecols=columns, chunksize=400_000):
                cgm = chunk.loc[chunk.PtID.isin(discovery.index) & chunk.RecordType.eq("CGM")].copy()
                if cgm.empty:
                    continue
                counts["selected_patient_cgm"] += len(cgm)
                utc = pd.to_datetime(cgm.UTCDtTm, format="ISO8601", utc=True)
                local = utc.dt.tz_convert(None) + cgm.PtID.map(offsets)
                counts["without_standard_hour"] += int(local.isna().sum())
                day = (local.dt.normalize() - cgm.PtID.map(dates)).dt.days
                slot = (60 * local.dt.hour + local.dt.minute) // 5
                in_window = day.between(1, landmark) | (day.between(landmark + 1, end) & slot.lt(72))
                if not in_window.any():
                    continue
                frame = cgm.loc[in_window, ["PtID", "CGMVal", "Units"]].copy()
                assert set(frame.Units.dropna()).issubset({"mmol/L", "mg/dL"})
                values = pd.to_numeric(frame.CGMVal, errors="raise") * frame.Units.map({"mmol/L": 18.01559, "mg/dL": 1})
                assert values.dropna().le(1000).all()
                counts["window_readings"] += len(frame)
                counts["invalid_or_missing"] += int((values.isna() | values.le(0)).sum())
                frame["glucose"] = values.where(values.gt(0)).clip(40, 400)
                counts["below_range"] += int((values.gt(0) & values.lt(40)).sum())
                counts["above_range"] += int(values.gt(400).sum())
                day = day.loc[in_window]
                frame["day"] = (day - landmark - day.le(landmark).astype(int)).astype("int16")
                frame["slot"] = slot.loc[in_window].astype("int16")
                frames.append(frame[["PtID", "day", "slot", "glucose"]].dropna())
        records = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=["PtID", "day", "slot", "glucose"])
        records.to_parquet(cache, index=False)
        metadata.write_text(json.dumps(counts, indent=2) + "\n")
        print(f"CGM{index} {counts['window_readings']}", flush=True)

records = pd.concat([pd.read_parquet(parts / f"device_{i}.parquet") for i in range(1, 7)], ignore_index=True)
assert set(records.PtID).issubset(discovery.index)
slots = records.groupby(["PtID", "day", "slot"], sort=True).glucose.median().reset_index()
assert slots.day.ne(0).all() and np.isfinite(slots.glucose).all()
assert slots.loc[slots.day.gt(0), "slot"].lt(72).all()
assert slots.day.min() >= -landmark and slots.day.max() <= end - landmark
slots.to_parquet(inputs / "discovery_slots.parquet", index=False)
counts = pd.DataFrame([json.loads((parts / f"device_{i}.json").read_text()) for i in range(1, 7)]).sum().to_dict()
counts.update({"valid_window_records": len(records), "unique_five_minute_slots": len(slots),
               "extra_records_in_occupied_slots": len(records) - len(slots), "reserved_patients_processed": 0,
               "clock": "UTC plus participant enrollment offset", "mmol_to_mg_dl": 18.01559})
(report / "source_counts.json").write_text(json.dumps(counts, indent=2) + "\n")
print(json.dumps(counts))
