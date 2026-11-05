import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.stats import norm

records = []
ratios = []
for tau in [0.25, 0.5, 0.75]:
    for scale in [0.5, 1.0, 2.0]:
        for row_scale in [0.5, 1.0, 2.0]:
            errors = []
            truth = scale * norm.ppf(tau)
            for bandwidth in [0.2, 0.1, 0.05, 0.025]:
                shift = (2 * tau - 1) * bandwidth / 2

                def equation(coefficient):
                    left = quad(lambda value: norm.cdf((coefficient + (value - shift) / row_scale) / scale),
                                -bandwidth, 0, epsabs=1e-12)[0]
                    right = quad(lambda value: norm.cdf((coefficient + (value - shift) / row_scale) / scale),
                                 0, bandwidth, epsabs=1e-12)[0]
                    return ((1 - tau) * left + tau * right) / bandwidth - tau

                target = brentq(equation, -20, 20, xtol=1e-13)
                mean = ((1 - tau) * (-bandwidth / 2 - shift) + tau * (bandwidth / 2 - shift)) / row_scale
                assert abs(mean) < 1e-12
                error = abs(target - truth)
                errors.append(error)
                records.append({"tau": tau, "outcome_scale": scale, "row_scale": row_scale,
                    "bandwidth": bandwidth, "quantile": truth, "smoothed_target": target,
                    "absolute_error": error, "kernel_mean": mean})
            if tau == 0.5:
                assert max(errors) < 1e-12
            else:
                ratio = errors[-2] / errors[-1]
                assert 3.8 < ratio < 4.2
                ratios.append(ratio)
output = Path("results/preflight/clinical_quantile_20261105")
pd.DataFrame(records).to_csv(output / "centering_population_checks.csv", index=False)
(output / "centering_acceptance.json").write_text(json.dumps({
    "status": "passed", "calculation": "deterministic_population_integration",
    "random_observations_generated": 0, "cases": len(records),
    "minimum_halving_error_ratio": min(ratios), "maximum_halving_error_ratio": max(ratios)}, indent=2) + "\n")
print(f"总体积分通过 {len(records)} 项检查，误差比范围 {min(ratios):.6f} 至 {max(ratios):.6f}")
