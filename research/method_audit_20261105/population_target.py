import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.stats import binom, norm

output = Path("results/preflight/method_audit_20261105")
output.mkdir(parents=True, exist_ok=True)


def cdf(value, components):
    return sum(weight * norm.cdf((value - mean) / sd) for weight, mean, sd in components)


def density(value, components):
    return sum(weight * norm.pdf((value - mean) / sd) / sd for weight, mean, sd in components)


def marginal_quantile(tau, components):
    return brentq(lambda value: cdf(value, components) - tau, -100, 100, xtol=1e-13)


def singleton_score(value, tau, penalty, components):
    score = 0.0
    lower, upper = (tau - 1) / penalty, tau / penalty
    for weight, mean, sd in components:
        left, right = (value + lower - mean) / sd, (value + upper - mean) / sd
        middle = norm.cdf(right) - norm.cdf(left)
        score += weight * (penalty * ((mean - value) * middle + sd * (norm.pdf(left) - norm.pdf(right)))
                           + (tau - 1) * norm.cdf(left) + tau * norm.sf(right))
    return score


def profile_score(value, tau, penalty, size, components):
    knots = (tau - 1 + np.arange(size + 1) / size) / penalty
    integral = 0.0
    for j in range(size):
        integral += quad(lambda shift: binom.sf(size - j - 1, size, cdf(value + shift, components)),
                         knots[j], knots[j + 1], epsabs=1e-11, epsrel=1e-11)[0]
    return penalty * (knots[-1] - integral)


def profile_target(tau, penalty, size, components):
    return brentq(lambda value: profile_score(value, tau, penalty, size, components),
                  -100, 100, xtol=1e-12)


def pinball_risk(value, tau, components):
    return sum(weight * (sd * norm.pdf((value - mean) / sd)
                         + (value - mean) * (norm.cdf((value - mean) / sd) - tau))
               for weight, mean, sd in components)


def singleton_objective(value, tau, penalty, components):
    lower, upper = (tau - 1) / penalty, tau / penalty

    def integrand(y):
        residual = y - value
        gamma = np.clip(residual, lower, upper)
        remainder = residual - gamma
        loss = remainder * (tau - float(remainder < 0)) + penalty * gamma**2 / 2
        return loss * density(y, components)

    points = [-np.inf, value + lower, value + upper, np.inf]
    return sum(quad(integrand, left, right, epsabs=1e-11, epsrel=1e-11)[0]
               for left, right in zip(points[:-1], points[1:]))


records = []
maximum_score_error = 0.0
standard = [(1.0, 0.0, 1.0)]
for size in [1, 2, 3, 6]:
    for tau in [0.25, 0.5, 0.75]:
        truth = marginal_quantile(tau, standard)
        for penalty in [0.5, 1, 2, 8]:
            target = profile_target(tau, penalty, size, standard)
            if size == 1:
                error = abs(profile_score(target, tau, penalty, size, standard)
                            - singleton_score(target, tau, penalty, standard))
                maximum_score_error = max(maximum_score_error, error)
            records.append({"cluster_size": size, "tau": tau, "Lambda": penalty,
                            "normal_quantile": truth, "profile_target": target,
                            "target_difference": target - truth,
                            "achieved_probability": cdf(target, standard),
                            "profile_prediction_loss": pinball_risk(target, tau, standard),
                            "quantile_prediction_loss": pinball_risk(truth, tau, standard)})
assert maximum_score_error < 1e-10
finite_difference_errors = []
for value in [-0.5, 0.2, 0.8]:
    step = 1e-5
    derivative = (singleton_objective(value + step, 0.75, 1, standard)
                  - singleton_objective(value - step, 0.75, 1, standard)) / (2 * step)
    finite_difference_errors.append(abs(derivative + singleton_score(value, 0.75, 1, standard)))
assert max(finite_difference_errors) < 1e-8
cases = {
    "conditional_scale": [(1, [(1.0, 0.0, 0.2)]), (1, [(1.0, 0.0, 2.0)])],
    "conditional_scale_repeated": [(2, [(1.0, 0.0, 0.2)]), (2, [(1.0, 0.0, 2.0)])],
    "cluster_size": [(1, standard), (3, standard)],
    "conditional_shape": [(1, [(0.8, 0.0, 0.05), (0.2, 4.0, 0.05)]), (1, standard)],
    "conditional_shape_repeated": [(2, [(0.8, 0.0, 0.05), (0.2, 4.0, 0.05)]), (2, standard)],
}
calibration_records = []
contrast_records = []
for case, groups in cases.items():
    tau, penalty = 0.75, 1.0
    targets = [profile_target(tau, penalty, size, components) for size, components in groups]
    truths = [marginal_quantile(tau, components) for size, components in groups]
    offset = brentq(lambda value: sum(cdf(target + value, components) for target, (_, components)
                                     in zip(targets, groups)) / 2 - tau, -20, 20, xtol=1e-12)
    probabilities = []
    for group, (target, truth, (size, components)) in enumerate(zip(targets, truths, groups)):
        probability = cdf(target + offset, components)
        probabilities.append(probability)
        calibration_records.append({"case": case, "group": group, "cluster_size": size, "tau": tau,
            "Lambda": penalty, "true_quantile": truth, "profile_target": target,
            "common_calibration": offset, "calibrated_prediction": target + offset,
            "conditional_probability": probability})
    assert abs(np.mean(probabilities) - tau) < 1e-10
    contrast_records.append({"case": case, "quantile_slope": truths[1] - truths[0],
                             "profile_slope": targets[1] - targets[0],
                             "calibrated_slope": targets[1] - targets[0],
                             "average_probability": float(np.mean(probabilities))})
pd.DataFrame(records).to_csv(output / "population_targets.csv", index=False)
pd.DataFrame(calibration_records).to_csv(output / "population_calibration.csv", index=False)
pd.DataFrame(contrast_records).to_csv(output / "population_slopes.csv", index=False)
checks = {"status": "passed", "calculation": "deterministic_population_integration",
          "random_samples_generated": 0, "clinical_data_used": False,
          "maximum_score_identity_error": maximum_score_error,
          "maximum_objective_derivative_error": max(finite_difference_errors)}
(output / "population_checks.json").write_text(json.dumps(checks, indent=2) + "\n")
print(pd.DataFrame(records).query("tau == 0.75 and Lambda == 1").to_string(index=False))
print(pd.DataFrame(contrast_records).to_string(index=False))
print(pd.DataFrame(calibration_records).to_string(index=False))
