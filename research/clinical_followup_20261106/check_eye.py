from pathlib import Path
import json
from zipfile import ZipFile

import numpy as np
import pandas as pd

from panel import study_config

config = study_config("uwhvf_central")
inputs = Path(config["input_directory"])
with ZipFile(Path(config["source_directory"]) / "release.zip") as archive:
    name, = [n for n in archive.namelist() if n.endswith("/alldata.json")]
    source = json.load(archive.open(name))
    name, = [n for n in archive.namelist() if n.endswith("/CSV/Coord_242.csv")]
    coords = pd.read_csv(archive.open(name))
central = np.flatnonzero((coords.X.abs() == 3) & (coords.Y.abs() == 3))
ring = np.flatnonzero(((coords.X.abs() == 3) & (coords.Y.abs() == 9)) |
                      ((coords.X.abs() == 9) & (coords.Y.abs() == 3)))
roles = pd.read_csv(inputs / "roles.csv", dtype={"patient": str}).set_index("patient")
rows = pd.read_csv(inputs / "discovery.csv").set_index("sample_id")
features = pd.read_csv(inputs / "features.csv.gz").set_index("sample_id")
audit = pd.read_csv(inputs / "time_audit.csv", dtype={"patient": str})
assert len(audit) == len(rows) == len(features) and audit.sample_id.is_unique
assert roles.loc[audit.patient, "role"].eq("discovery").all()
errors = []
target_errors = []


def field_at(series, age, key):
    visits = [v for v in series if abs(v["age"] - age) <= 0.00006]
    assert visits, age
    return np.median(np.array([v[key] for v in visits], dtype=float), axis=0)


for record in audit.itertuples():
    series = source["data"][record.patient][record.eye[0]]
    assert record.first_age < record.landmark < record.outcome_age
    future = field_at(series, record.outcome_age, "hvf_seq")[central]
    errors.append(abs(future.mean() - rows.loc[record.sample_id, "y"]))
    td = np.stack([field_at(series, t, "td_seq")[ring] for t in [record.first_age, record.landmark]])
    original = features.loc[record.sample_id, "C_ring_persistent"]
    if np.isfinite(original):
        # CSV 公开的总差值保留两位小数。
        value = (np.round(td, 2) <= -5).all(axis=0).sum() / 2
        target_errors.append(abs(value - original))
maximum = float(max(errors))
assert maximum <= 0.006 and max(target_errors) == 0, (maximum, max(target_errors))
state = {"source": "independent published JSON arrays", "patients": audit.patient.nunique(), "records": len(audit),
         "clinical_coefficients_checked_rows": len(target_errors), "maximum_endpoint_difference_db": maximum,
         "maximum_target_input_difference": float(max(target_errors)), "reserved_outcomes_used": False}
(Path(config["report_directory"]) / "endpoint_verification.json").write_text(json.dumps(state, indent=2) + "\n")
print(json.dumps(state))
