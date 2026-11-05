import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

config = yaml.safe_load(Path("config/application_charite_hb_20261105.yml").read_text())
root = Path(config["output_root"])
report = Path(config["report_directory"])
training = pd.read_csv(Path(config["input_directory"]) / "training.csv")
validation = pd.read_csv(Path(config["input_directory"]) / "validation.csv")
assert not set(training.strain_id) & set(validation.strain_id)
assert training.groupby("strain_id").inner_fold.nunique().eq(1).all()
prediction = pd.concat([pd.read_csv(root / "validation" / "predictions_R.csv"),
                        pd.read_csv(root / "validation" / "predictions_boost.csv")])
prediction = prediction.merge(validation[["sample_id", "strain_id", "y"]],
                              on=["sample_id", "strain_id"], validate="many_to_one")
assert len(prediction) == len(validation) * 3 * 11
wide = prediction.pivot(index=["method", "sample_id", "strain_id", "y"], columns="tau", values="prediction")
assert list(wide.columns) == config["quantiles"] and wide.notna().all().all()
ordered = wide.copy()
ordered.iloc[:, :] = np.sort(np.maximum(wide.to_numpy(), 0), axis=1)
ordered = ordered.reset_index().melt(id_vars=["method", "sample_id", "strain_id", "y"],
                                    var_name="tau", value_name="prediction")
records = []
for scoring, frame in [("raw", prediction), ("primary", ordered)]:
    frame = frame.copy()
    residual = frame.y - frame.prediction
    frame["loss"] = np.where(residual >= 0, frame.tau * residual, (frame.tau - 1) * residual)
    frame["scoring"] = scoring
    records.append(frame)
reconstructed = pd.concat(records, ignore_index=True)
saved = pd.read_csv(root / "validation" / "scored_predictions.csv")
keys = ["scoring", "method", "sample_id", "strain_id", "tau"]
comparison = reconstructed.merge(saved, on=keys, validate="one_to_one", suffixes=("_new", "_saved"))
assert len(comparison) == len(saved) == 2 * len(prediction)
errors = {field: float(np.max(np.abs(comparison[f"{field}_new"] - comparison[f"{field}_saved"])))
          for field in ["y", "prediction", "loss"]}
assert max(errors.values()) < 1e-10
patient_loss = reconstructed.groupby(["scoring", "method", "tau", "strain_id"], observed=True).loss.mean()
summary = patient_loss.groupby(["scoring", "method", "tau"], observed=True).agg(["mean", "std", "count"]).reset_index()
summary["se"] = summary["std"] / np.sqrt(summary["count"])
saved_summary = pd.read_csv(report / "prediction_summary.csv")
summary = summary.merge(saved_summary, on=["scoring", "method", "tau"], validate="one_to_one", suffixes=("_new", "_saved"))
errors["patient_loss"] = float((summary["mean"] - summary.loss).abs().max())
errors["patient_standard_error"] = float((summary.se_new - summary.se_saved).abs().max())
assert summary["count"].eq(validation.strain_id.nunique()).all()
comparisons = pd.read_csv(report / "prediction_comparisons.csv")
gains = []
for row in comparisons.itertuples():
    reference = patient_loss.loc[(row.scoring, row.comparator, row.tau)]
    target = patient_loss.loc[(row.scoring, row.target_method, row.tau)]
    assert reference.index.equals(target.index)
    absolute = float((reference - target).mean())
    relative = absolute / float(reference.mean())
    gains.append(max(abs(absolute - row.absolute_gain), abs(relative - row.relative_gain)))
errors["paired_gain"] = max(gains)
assert max(errors.values()) < 1e-10
assert not comparisons.query("scoring == 'primary' and tau == 0.25").superiority_pass.any()
release = json.loads((root / "validation_release.json").read_text())
assert len(release["model_directories"]) == 19
status = {"status": "passed", "training_patients": training.strain_id.nunique(),
          "validation_patients": validation.strain_id.nunique(), "validation_records": len(validation),
          "methods": prediction.method.nunique(), "prediction_rows_checked": len(saved),
          "maximum_errors": errors, "patient_overlap": 0}
(report / "python_audit_status.json").write_text(json.dumps(status, indent=2) + "\n")
print(json.dumps(status, indent=2))
