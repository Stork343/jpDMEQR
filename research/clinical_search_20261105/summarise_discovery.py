import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

parser = argparse.ArgumentParser()
parser.add_argument("config", type=Path)
parser.add_argument("--comparators-only", action="store_true")
args = parser.parse_args()
config = yaml.safe_load(args.config.read_text())
root = Path(config["output_root"])
output = Path(config["report_directory"])
progress = pd.read_csv(root / "progress.csv")
assert len(progress) == config["outer_folds"] * 15
assert progress.status.eq("completed").all() and progress.exit_code.dropna().eq(0).all()
rows = pd.read_csv(Path(config["input_directory"]) / "discovery.csv")
frames = []
runtime = []
for job in progress.itertuples():
    directory = Path(job.directory)
    status = json.loads((directory / "status.json").read_text())
    assert status["status"] == "completed" and not status["formal_validation_used"]
    prediction = pd.read_csv(directory / "prediction.csv")
    expected = rows.loc[rows.outer_fold.eq(int(job.trait.removeprefix("outer_")))]
    assert set(prediction.sample_id) == set(expected.sample_id)
    frames.append(prediction)
    runtime.append({"outer_fold": job.trait, "method": job.method, "runtime_sec": status["runtime_sec"]})
prediction = pd.concat(frames, ignore_index=True)
keys = ["sample_id", "strain_id", "method", "tau"]
assert not prediction[keys].duplicated().any()
assert len(prediction) == len(rows) * 7 * 3
prediction = prediction.merge(rows[["sample_id", "strain_id", "y", "outer_fold"]],
                              on=["sample_id", "strain_id"], validate="many_to_one")
wide = prediction.pivot(index=["method", "sample_id", "strain_id", "y", "outer_fold"], columns="tau", values="prediction")
assert list(wide.columns) == config["quantiles"] and wide.notna().all().all()
ordered = wide.copy()
ordered.iloc[:, :] = np.sort(np.maximum(wide.to_numpy(), 0), axis=1)
ordered = ordered.reset_index().melt(id_vars=["method", "sample_id", "strain_id", "y", "outer_fold"],
                                    var_name="tau", value_name="prediction")
scored = []
for scoring, current in [("raw", prediction), ("primary", ordered)]:
    current = current.copy()
    current["scoring"] = scoring
    error = current.y - current.prediction
    current["loss"] = np.maximum(current.tau * error, (current.tau - 1) * error)
    scored.append(current)
scored = pd.concat(scored, ignore_index=True)
patient = scored.groupby(["scoring", "method", "tau", "strain_id", "outer_fold"], observed=True).loss.mean().reset_index()
summary = patient.groupby(["scoring", "method", "tau"], observed=True).agg(
    patients=("loss", "size"), loss=("loss", "mean"), sd=("loss", "std")).reset_index()
summary["se"] = summary.sd / np.sqrt(summary.patients)
summary.to_csv(output / "discovery_loss.csv", index=False)
folds = patient.groupby(["scoring", "method", "tau", "outer_fold"], observed=True).loss.mean().reset_index()
folds.to_csv(output / "discovery_fold_loss.csv", index=False)
pd.DataFrame(runtime).to_csv(output / "discovery_runtime.csv", index=False)
scored.to_csv(root / "discovery_scored_predictions.csv", index=False)
if args.comparators_only:
    print(summary.loc[summary.scoring.eq("primary")].pivot(index="method", columns="tau", values="loss").to_string())
    raise SystemExit(0)
comparators = [method for method in config["comparators"] if not method.startswith(("QRF", "XGB", "LGB"))]
records = []
for scoring in ["primary", "raw"]:
    for tau in config["quantiles"]:
        risks = patient.loc[patient.scoring.eq(scoring) & patient.tau.eq(tau)].pivot(
            index="strain_id", columns="method", values="loss")
        assert risks.notna().all().all() and len(risks) == rows.strain_id.nunique()
        for target in config["target_candidates"]:
            for comparator in comparators:
                reference = float(risks[comparator].mean())
                proposed = float(risks[target].mean())
                records.append({"scoring": scoring, "tau": tau, "target": target, "comparator": comparator,
                                "target_loss": proposed, "comparator_loss": reference,
                                "relative_gain": (reference - proposed) / reference,
                                "formal_validation_used": False})
comparison = pd.DataFrame(records)
comparison.to_csv(output / "discovery_comparisons.csv", index=False)
decision = comparison.loc[comparison.scoring.eq("primary") & comparison.tau.eq(config["decision_quantile"])]
minimum = decision.groupby("target").relative_gain.min()
selected = str(minimum.idxmax())
passed = bool(minimum.loc[selected] >= config["screening_minimum_relative_gain"])
status = {"study_id": config["study_id"], "status": "linear_screen_completed", "completed_jobs": len(progress),
          "observed_process_exit_codes": int(progress.exit_code.notna().sum()),
          "completed_artifacts_without_exit_code": int(progress.exit_code.isna().sum()),
          "discovery_patients": rows.strain_id.nunique(), "discovery_records": len(rows),
          "decision_quantile": config["decision_quantile"], "best_screening_candidate": selected,
          "minimum_gain_by_candidate": minimum.to_dict(), "continue_nonlinear_comparison": passed,
          "formal_validation_used": False, "manuscript_adoption": False}
(output / "discovery_status.json").write_text(json.dumps(status, indent=2) + "\n")
print(summary.loc[summary.scoring.eq("primary")].pivot(index="method", columns="tau", values="loss").to_string())
print(json.dumps(status, indent=2))
