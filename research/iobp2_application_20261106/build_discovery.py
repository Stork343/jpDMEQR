import json
from pathlib import Path
import sys
from zipfile import ZipFile

import numpy as np
import pandas as pd
import yaml

config = yaml.safe_load(Path(sys.argv[1] if len(sys.argv) > 1 else "config/application_iobp2_20261106.yml").read_text())
inputs = Path(config["input_directory"])
report = Path(config["report_directory"])
roles = pd.read_csv(inputs / "roles.csv", dtype={"PtID": str})
discovery = roles.loc[roles.role.eq("discovery")].set_index("PtID")
dates = pd.to_datetime(discovery.RandDt.str.split().str[0], format="%m/%d/%Y")
slots = pd.read_parquet(inputs / "discovery_slots.parquet")
assert set(slots.PtID).issubset(discovery.index)
assert not slots.duplicated(["PtID", "day", "slot"]).any()


def tail_mean(series):
    values = np.sort(series.to_numpy())[::-1]
    size, remainder = divmod(len(values), 10)
    return (values[:size].sum() + remainder / 10 * values[size]) / (len(values) / 10)


baseline = slots.loc[slots.day.lt(0)]
if config.get("baseline_window") == "run_in":
    pre_days = baseline.groupby(["PtID", "day"]).size().rename("n").reset_index()
    latest = (pre_days.loc[pre_days.n.ge(config["minimum_slots_per_day"])]
              .sort_values(["PtID", "day"]).groupby("PtID").tail(config["minimum_baseline_days"]))
    baseline = baseline.merge(latest[["PtID", "day"]], validate="many_to_one")
    slots = pd.concat([baseline, slots.loc[slots.day.gt(0)]], ignore_index=True)
    lag = latest.groupby("PtID").day.agg(first="min", last="max", days="size")
    lag.loc[lag.days.eq(config["minimum_baseline_days"]), ["first", "last"]].mul(-1).describe().to_csv(
        report / "baseline_recency_days.csv")
night_slots = slots.loc[slots.slot.lt(72)]
nightly = night_slots.groupby(["PtID", "day"]).glucose.agg(count="count", tail=tail_mean)
nightly = nightly.loc[nightly["count"].ge(config["minimum_slots_per_night"])]
baseline_nights = nightly.loc[nightly.index.get_level_values("day") < 0]
days = baseline.groupby(["PtID", "day"]).size()
qualified_days = days.loc[days.ge(config["minimum_slots_per_day"])].groupby(level="PtID").size()
coverage = pd.DataFrame(index=discovery.index)
coverage["baseline_days"] = qualified_days
coverage["baseline_nights"] = baseline_nights.groupby(level="PtID").size()
future = nightly.reset_index().query("day > 0").copy()
future["period"] = (future.day - 1) // config["period_days"] + 1
periods = future.groupby(["PtID", "period"])["tail"]
for period in range(1, config["periods"] + 1):
    coverage[f"period_{period}_nights"] = periods.size().unstack().get(period)
coverage = coverage.fillna(0).astype(int)
baseline_ok = ((coverage.baseline_days >= config["minimum_baseline_days"])
               & (coverage.baseline_nights >= config["minimum_baseline_nights"]))
all_periods = coverage.filter(regex="^period_").ge(config["minimum_nights_per_period"]).all(axis=1)
eligible = coverage.index[baseline_ok & all_periods]
inventory = coverage.join(discovery[["TrtGroup", "outer_fold"]])
inventory["baseline_eligible"] = baseline_ok
inventory["three_periods_eligible"] = all_periods
inventory["analysis_eligible"] = baseline_ok & all_periods
inventory.to_csv(inputs / "discovery_coverage.csv")
by_arm = inventory.groupby("TrtGroup").agg(randomized=("outer_fold", "size"),
    baseline_eligible=("baseline_eligible", "sum"), three_periods_eligible=("three_periods_eligible", "sum"),
    analysis_eligible=("analysis_eligible", "sum"))
by_arm.to_csv(report / "coverage_by_arm.csv")
by_fold = inventory.groupby("outer_fold").agg(randomized=("TrtGroup", "size"),
    analysis_eligible=("analysis_eligible", "sum"))
by_fold.to_csv(report / "coverage_by_fold.csv")
outcomes = periods.mean().loc[lambda x: x.index.get_level_values("PtID").isin(eligible)]
assert len(outcomes) == len(eligible) * config["periods"] and outcomes.between(40, 400).all()
boundary = float(outcomes.eq(400).mean()) if len(outcomes) else None
checks = {"patients": bool(len(eligible) >= config["minimum_discovery_patients"]),
          "arms": bool(by_arm.analysis_eligible.min() >= config["minimum_per_arm"]),
          "folds": bool(by_fold.analysis_eligible.min() >= config["minimum_per_fold"]),
          "endpoint_support": outcomes.nunique() >= 20, "boundary": boundary is not None and boundary < 0.05}
status = {"discovery_patients": len(eligible), "period_records": len(outcomes),
          "distinct_endpoint_values": outcomes.nunique(), "periods_at_sensor_ceiling": boundary,
          "criteria": checks, "passed": all(checks.values()), "formal_validation_used": False,
          "manuscript_adoption": False}
(report / "eligibility.json").write_text(json.dumps(status, indent=2) + "\n")
print(json.dumps(status), flush=True)
if not status["passed"]:
    raise SystemExit(0)

base = discovery.loc[eligible]
features = pd.DataFrame(index=eligible)
features["C_age"] = base.AgeAsofEnrollDt
features["C_pediatric"] = base.pediatric.astype(int)
features["C_armBP"] = base.TrtGroup.eq("BP").astype(int)

with ZipFile(Path(config["source_directory"]) / "IOBP2_RCT.zip") as archive:
    def read_fields(table, fields):
        with archive.open(f"Data Tables/{table}.txt") as stream:
            frame = pd.read_csv(stream, sep="|", dtype=str, usecols=["PtID", *fields])
        return frame.loc[frame.PtID.isin(eligible)].copy()

    sex = read_fields("IOBP2DiabScreening", ["Sex"]).set_index("PtID").Sex
    assert sex.index.is_unique and set(sex.dropna()).issubset({"M", "F"})
    features["C_female"] = sex.map({"F": 1, "M": 0})
    modes = read_fields("IOBP2RandBaseInfo", ["DiabTrt"]).set_index("PtID").DiabTrt.reindex(eligible)
    assert modes.index.is_unique
    categories = {"Conventional pump without automation": "pump", "Patch pump": "patch",
                  "Tandem Control IQ": "control_iq", "Tandem Basal IQ": "basal_iq",
                  "Medtronic 670G using automode": "medtronic", "Medtronic 670G/770G using automode": "medtronic",
                  "Multiple daily injections": "injections"}
    assert set(modes.dropna()).issubset(categories)
    for code in sorted(set(categories.values()) - {"injections"}):
        features[f"C_mode_{code}"] = modes.map(categories).eq(code).astype(float).where(modes.notna())
    local = read_fields("IOBP2DiabLocalHbA1c", ["Visit", "HbA1cTestDt", "HbA1cTestRes"])
    local["available"] = pd.to_datetime(local.HbA1cTestDt, format="%m/%d/%Y")
    local = local.loc[local.Visit.eq("Screening") & local.available.le(local.PtID.map(dates))
                      & local.HbA1cTestRes.notna(), ["PtID", "available", "HbA1cTestRes"]]
    central = read_fields("STASampleResults", ["Visit", "CollectionDt", "AnalysisDt", "ResultName", "Value",
                                              "Units", "STAResultStatus"])
    central = central.loc[central.Visit.eq("Randomization") & central.ResultName.eq("GLYHB")
                          & central.Units.eq("%") & central.STAResultStatus.eq("Complete")].copy()
    central["available"] = pd.to_datetime(central.AnalysisDt, format="%m/%d/%Y")
    collection = pd.to_datetime(central.CollectionDt, format="%m/%d/%Y")
    # 同日仅有日期时，中央检验报告须在随机日前已有结果。
    central = central.loc[collection.le(central.available) & central.available.lt(central.PtID.map(dates)),
                          ["PtID", "available", "Value"]].rename(columns={"Value": "HbA1cTestRes"})
    known = pd.concat([local, central]).sort_values("available").drop_duplicates("PtID", keep="last")
    features["C_hba1c"] = pd.to_numeric(known.set_index("PtID").HbA1cTestRes, errors="raise")

b = baseline.loc[baseline.PtID.isin(eligible)].copy()
g = b.groupby("PtID").glucose
features["C_meanCGM"] = g.mean()
features["C_cv"] = 100 * g.std() / g.mean()
features["C_TIR70_180"] = 100 * b.glucose.between(70, 180).groupby(b.PtID).mean()
features["C_TBR70"] = 100 * b.glucose.lt(70).groupby(b.PtID).mean()
features["C_TAR250"] = 100 * b.glucose.gt(250).groupby(b.PtID).mean()
features["C_baseline_night_tail"] = baseline_nights.groupby(level="PtID")["tail"].mean()
bn = b.loc[b.slot.lt(72)].merge(baseline_nights.reset_index()[["PtID", "day"]], on=["PtID", "day"])
half = bn.groupby(["PtID", "day", bn.slot.ge(36)]).glucose.mean().unstack()
assert half.notna().all().all()
features["C_dawn"] = (half[True] - half[False]).groupby(level="PtID").mean() / 20
features["C_arm_dawn"] = features.C_dawn * features.C_armBP
b["bin"] = b.slot // 3
by_bin = b.groupby(["PtID", "bin"])
profile = by_bin.glucose.quantile([0.1, 0.5, 0.9]).unstack()
profile.loc[by_bin.day.nunique() < config["minimum_days_per_profile_bin"]] = np.nan
profile = profile.unstack("bin").reindex(columns=pd.MultiIndex.from_product([[0.1, 0.5, 0.9], range(96)]))
profile.columns = [f"P_q{int(q * 100):02d}_bin{k:02d}" for q, k in profile.columns]
features = features.join(profile)
assert profile.shape[1] == 288 and features.C_arm_dawn.notna().all()
rows = outcomes.rename("y").reset_index()
rows["sample_id"] = "IO_" + rows.PtID + "_P" + rows.period.astype(str)
rows["strain_id"] = "IO_" + rows.PtID
rows["outer_fold"] = rows.PtID.map(base.outer_fold)
rows["condition"] = config["response"]
rows["role"] = "discovery"
matrix = rows[["sample_id", "PtID", "period"]].join(features, on="PtID", validate="many_to_one")
matrix["C_period2"] = matrix.period.eq(2).astype(int)
matrix["C_period3"] = matrix.period.eq(3).astype(int)
assert matrix.sample_id.is_unique and rows.groupby("PtID").outer_fold.nunique().eq(1).all()
matrix.drop(columns=["PtID", "period"]).to_csv(inputs / "features.csv.gz", index=False)
rows[["sample_id", "strain_id", "y", "condition", "role", "outer_fold"]].to_csv(inputs / "discovery.csv", index=False)
