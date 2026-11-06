import json
from pathlib import Path
import sys

import duckdb
import numpy as np
import pandas as pd
import yaml

config = yaml.safe_load(Path(sys.argv[1]).read_text())
inputs = Path(config["input_directory"])
report = Path(config["report_directory"])
status = json.loads((report / "eligibility.json").read_text())
connection = duckdb.connect()
connection.read_parquet(str(inputs / "discovery_slots.parquet")).create_view("slots")
roles = pd.read_csv(inputs / "roles.csv", dtype={"PtID": str})
assert roles.PtID.is_unique
if "SiteID" in roles:
    assert roles.groupby("SiteID").role.nunique().eq(1).all()
    assert roles.groupby("SiteID").outer_fold.nunique().eq(1).all()
actual_ids = set(connection.sql("select distinct PtID from slots").df().PtID)
assert actual_ids.issubset(set(roles.loc[roles.role.eq("discovery"), "PtID"]))
days_needed = config["minimum_baseline_days"]
occupied = config["minimum_slots_per_day"]
baseline = f"""
    select PtID, day from slots where day < 0 group by PtID, day having count(*) >= {occupied}
"""
if config.get("baseline_window") == "run_in":
    baseline += f" qualify row_number() over (partition by PtID order by day desc) <= {days_needed}"
connection.sql(baseline).create_view("eligible_days")
connection.sql("""
    select * from slots where day > 0 or exists
      (select 1 from eligible_days d where d.PtID=slots.PtID and d.day=slots.day)
""").create_view("selected")
night_source = "selected" if config.get("baseline_window") == "run_in" else "slots"
connection.sql(f"""
    with ranked as (
      select *, row_number() over (partition by PtID,day order by glucose desc, slot) as position,
             count(*) over (partition by PtID,day) as n from {night_source} where slot < 72
    )
    select PtID, day, sum(glucose * least(1, greatest(0, n * 0.1 - position + 1))) / (0.1*n) as tail,
           avg(glucose) filter (where slot >= 36) - avg(glucose) filter (where slot < 36) as dawn
    from ranked where n >= {config["minimum_slots_per_night"]} group by PtID,day,n
""").create_view("nights")
stats = connection.sql(f"""
    with p as (select PtID, floor((day - 1)/{config["period_days"]})::integer+1 as period,
                        count(*) as n, avg(tail) as y
               from nights where day > 0 group by PtID,period),
    satisfying as (select PtID from p group by PtID having count(*)={config["periods"]}
                   and min(n)>={config["minimum_nights_per_period"]}),
    b as (select PtID, avg(tail) as night_tail, avg(dawn)/20 as dawn from nights where day < 0
          group by PtID having count(*) >= {config["minimum_baseline_nights"]}),
    d as (select PtID from eligible_days group by PtID having count(*) >= {days_needed})
    select p.PtID, period, y, night_tail, dawn from p join satisfying using(PtID) join b using(PtID)
      join d using(PtID) order by p.PtID, period
""").df()
assert stats.PtID.nunique() == status["discovery_patients"] and len(stats) == status["period_records"]
coverage = pd.read_csv(inputs / "discovery_coverage.csv", dtype={"PtID": str})
assert set(stats.PtID) == set(coverage.loc[coverage.analysis_eligible, "PtID"])
audit = {"discovery_patients_reconstructed": stats.PtID.nunique(), "records_reconstructed": len(stats),
         "reserved_patient_readings": len(actual_ids & set(roles.loc[roles.role.eq("reserved"), "PtID"]))}
if status["passed"]:
    rows = pd.read_csv(inputs / "discovery.csv")
    features = pd.read_csv(inputs / "features.csv.gz")
    stats["sample_id"] = config.get("patient_prefix", "IO") + "_" + stats.PtID + "_P" + stats.period.astype(str)
    paired = stats.merge(rows, on="sample_id", validate="one_to_one", suffixes=("_SQL", "_Python"))
    paired = paired.merge(features, on="sample_id", validate="one_to_one")
    assert len(paired) == len(rows) == len(stats) == len(features)
    audit["maximum_endpoint_difference"] = float(abs(paired.y_SQL - paired.y_Python).max())
    audit["maximum_baseline_night_tail_difference"] = float(abs(paired.night_tail - paired.C_baseline_night_tail).max())
    audit["maximum_dawn_difference"] = float(abs(paired.dawn - paired.C_dawn).max())
    assert max(audit[key] for key in audit if key.startswith("maximum")) < 1e-9
    profiles = connection.sql(f"""
        select PtID, (slot // 3) as bin,
               quantile_cont(glucose, [0.1,0.5,0.9]) as q from {night_source} where day < 0
          group by PtID,bin having count(distinct day) >= {config["minimum_days_per_profile_bin"]}
    """).df().explode("q")
    profiles["quantile"] = profiles.groupby(["PtID", "bin"]).cumcount().map({0: 10, 1: 50, 2: 90})
    profiles["name"] = [f"P_q{q:02d}_bin{k:02d}" for q, k in zip(profiles["quantile"], profiles.bin)]
    check = paired.loc[paired.period.eq(1)].set_index("PtID").filter(regex="^P_")
    reference = profiles.pivot(index="PtID", columns="name", values="q").reindex_like(check).astype(float)
    np.testing.assert_allclose(reference, check, rtol=1e-12, atol=1e-10, equal_nan=True)
    assert paired.groupby("PtID").outer_fold.nunique().eq(1).all()
    assert paired.groupby("PtID")[list(check.columns)].nunique(dropna=False).eq(1).all().all()
    audit["period_prediction_rows"] = len(rows)
    audit["profile_columns_reconstructed"] = len(check.columns)
(report / "panel_verification.json").write_text(json.dumps(audit, indent=2) + "\n")
print(json.dumps(audit))
