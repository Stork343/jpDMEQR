import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

config = yaml.safe_load(Path("config/clinical_quantile_20261105.yml").read_text())
output = Path("results/preflight/clinical_quantile_20261105")
summaries = []
comparisons = []
checks = []
for path in config["source_studies"]:
    study = yaml.safe_load(Path(path).read_text())
    root = Path(config["output_root"]) / study["study_id"]
    status = json.loads((root / "status.json").read_text())
    if status["status"] != "completed":
        continue
    predictions = pd.read_csv(root / "predictions.csv")
    assert not predictions.duplicated(["sample_id", "tau"]).any()
    assert predictions.groupby("strain_id")["fold"].nunique().eq(1).all()
    wide = predictions.pivot(index="sample_id", columns="tau", values="prediction")
    assert list(wide.columns) == config["quantiles"] and not wide.isna().any().any()
    monotone = pd.DataFrame(np.maximum(0, np.sort(wide.to_numpy(), axis=1)),
                           index=wide.index, columns=wide.columns)
    monotone = monotone.rename_axis(columns="tau").stack().rename("sorted_prediction").reset_index()
    predictions = predictions.merge(monotone, on=["sample_id", "tau"], validate="one_to_one")
    for rule, column in [("raw", "prediction"), ("nonnegative_sorted_quantiles", "sorted_prediction")]:
        residual = predictions["y"] - predictions[column]
        predictions["loss"] = residual * (predictions["tau"] - (residual < 0))
        patient = predictions.groupby(["tau", "strain_id"])["loss"].mean().reset_index()
        summary = patient.groupby("tau")["loss"].agg(["mean", "count"]).reset_index()
        summary["study_id"] = study["study_id"]
        summary["method"] = "CQ_ENET"
        summary["rule"] = rule
        summaries.append(summary.rename(columns={"mean": "loss", "count": "patients"}))
    archived_path = Path(study["report_directory"]) / "discovery_loss.csv"
    archived = pd.read_csv(archived_path)
    assert {"method", "tau", "loss"}.issubset(archived.columns)
    archived = archived[archived["scoring"].eq("primary")]
    assert not archived.duplicated(["method", "tau"]).any()
    current = summary.set_index("tau")["mean"]
    for row in archived.itertuples():
        comparisons.append({"study_id": study["study_id"], "tau": row.tau,
            "method": "CQ_ENET", "comparator": row.method, "current_loss": current.loc[row.tau],
            "comparator_loss": row.loss, "relative_improvement": 1 - current.loc[row.tau] / row.loss,
            "decision_quantile": row.tau == study["decision_quantile"]})
    diagnostics = pd.read_csv(root / "diagnostics.csv")
    assert len(diagnostics) == study["outer_folds"] * len(config["quantiles"])
    assert diagnostics["gap"].abs().max() < config["exact_fit_tolerance"]
    assert diagnostics["path_objective_difference"].abs().max() < config["maximum_path_objective_difference"]
    diagnostics.to_csv(output / f"{study['study_id']}_fit_diagnostics.csv", index=False)
    checks.append({"study_id": study["study_id"], "fits": len(diagnostics),
                  "maximum_gap": diagnostics["gap"].abs().max(),
                  "maximum_stationarity": diagnostics["stationarity"].max(),
                  "total_fitting_seconds": diagnostics["runtime_sec"].sum()})
assert summaries, "尚无完整队列结果"
pd.concat(summaries, ignore_index=True).to_csv(output / "discovery_losses.csv", index=False)
pd.DataFrame(comparisons).to_csv(output / "discovery_comparisons.csv", index=False)
pd.DataFrame(checks).to_csv(output / "discovery_numerical_checks.csv", index=False)
(output / "discovery_status.json").write_text(json.dumps({
    "completed_studies": [record["study_id"] for record in checks],
    "formal_validation_used": False, "manuscript_adoption": False}, indent=2) + "\n")
print(pd.DataFrame(comparisons).query("decision_quantile").to_string(index=False))
