import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.stats import norm
import yaml

config = yaml.safe_load(Path(sys.argv[1] if len(sys.argv) > 1 else "config/clinical_quantile_20261105.yml").read_text())
output = Path(config["report_directory"])
records = []
for directory in sorted(Path(config["output_root"]).glob("orthogonal_*")):
    status = json.loads((directory / "status.json").read_text())
    assert status["status"] == "completed"
    patients = pd.read_csv(directory / "patient_influence.csv")
    assert not patients.duplicated(["method", "strain_id"]).any()
    wide = patients.pivot(index=["strain_id", "fold"], columns="method", values="centered")
    assert not wide.isna().any().any()
    values = wide[["independence", "exchangeable"]].to_numpy()
    covariance = values.T @ values / len(values)
    difference_variance = np.mean((values[:, 1] - values[:, 0]) ** 2)
    assert difference_variance >= 0
    weight = 0.0 if difference_variance == 0 else float(np.clip(
        (covariance[0, 0] - covariance[0, 1]) / difference_variance, 0, 1))
    weights = np.array([1 - weight, weight])
    points = patients.assign(corrected=patients["initial"] + patients["influence"]).groupby("method")["corrected"].mean()
    estimate = weights @ points.loc[["independence", "exchangeable"]].to_numpy()
    influence = values @ weights
    variance = np.mean(influence**2)
    assert variance <= min(np.diag(covariance)) + 1e-12
    grid = np.linspace(0, 1, 1001)
    grid_variance = (1 - grid)**2 * covariance[0, 0] + grid**2 * covariance[1, 1] + 2 * grid * (1 - grid) * covariance[0, 1]
    assert variance <= grid_variance.min() + 1e-12
    standard_error = np.sqrt(variance / len(values))
    squared = influence**2
    assert squared.sum() > 0
    contributions = np.sort(squared / squared.sum())
    record = {"study_id": directory.name.removeprefix("orthogonal_"), "patients": len(values),
              "target": status["target"], "tau": status["tau"], "exchangeable_weight": weight,
              "estimate": float(estimate), "standard_error": standard_error,
              "lower": estimate - norm.ppf(0.975) * standard_error,
              "upper": estimate + norm.ppf(0.975) * standard_error,
              "independence_standard_error": np.sqrt(covariance[0, 0] / len(values)),
              "exchangeable_standard_error": np.sqrt(covariance[1, 1] / len(values)),
              "estimated_se_reduction": 1 - standard_error / np.sqrt(covariance[0, 0] / len(values)),
              "largest_patient_variance_share": float(contributions[-1]),
              "largest_five_patient_variance_share": float(contributions[-5:].sum()),
              "formal_validation_used": False, "manuscript_adoption": False}
    records.append(record)
assert records
pd.DataFrame(records).to_csv(output / "adaptive_inference.csv", index=False)
(output / "adaptive_inference_checks.json").write_text(json.dumps({
    "status": "passed", "checks": ["patient_pairing", "covariance_identity", "convex_weight_minimisation"],
    "role": "exploratory_estimator_development", "formal_validation_used": False}, indent=2) + "\n")
print(pd.DataFrame(records).to_string(index=False))
