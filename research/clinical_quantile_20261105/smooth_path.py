import time

import numpy as np
from skglm.experimental import QuantileHuber
from skglm.penalties import WeightedL1_plus_L2
from skglm.solvers import AndersonCD


def fit_smooth_path(X, y, row_weights, penalty_factors, lambdas, alpha, tau, patients, initial_beta, tolerance=1e-6):
    factors = np.asarray(penalty_factors, dtype=np.float64).reshape(-1).copy()
    weights = len(y) * np.asarray(row_weights, dtype=np.float64).reshape(-1) / patients
    original_X = np.asfortranarray(np.asarray(X, dtype=np.float64) * weights[:, None])
    bandwidth = float(patients)**(-0.2)
    response = np.asarray(y, dtype=np.float64).reshape(-1) * weights + (2 * tau - 1) * bandwidth / 2
    clinical = np.flatnonzero(factors == 0)
    penalized = np.flatnonzero(factors > 0)
    Q, R = np.linalg.qr(original_X[:, clinical], mode="reduced")
    assert np.min(np.abs(np.diag(R))) > 1e-10
    projection = np.linalg.solve(R, Q.T @ original_X[:, penalized])
    transformed = original_X.copy(order="F")
    transformed[:, clinical] = Q * np.sqrt(len(y))
    transformed[:, penalized] -= original_X[:, clinical] @ projection
    retained = (factors == 0) | (np.linalg.norm(transformed, axis=0) > 1e-10)
    design = np.asfortranarray(transformed[:, retained])
    initial = np.asarray(initial_beta, dtype=np.float64).reshape(-1).copy()
    initial[clinical] = R @ (initial[clinical] + projection @ initial[penalized]) / np.sqrt(len(y))
    coefficients = initial[retained].copy()
    numerical_ridge = 1 / float(patients)
    datafit = QuantileHuber(quantile=float(tau), delta=bandwidth)
    solver = AndersonCD(max_iter=100, max_epochs=20000, p0=max(10, 2 * len(clinical)), tol=float(tolerance),
                        fit_intercept=False, warm_start=True, verbose=0)
    beta = []
    diagnostics = []
    for index, strength in enumerate(np.asarray(lambdas).reshape(-1)):
        started = time.perf_counter()
        penalty = WeightedL1_plus_L2(1.0, 0.5,
            2 * float(strength) * float(alpha) * factors[retained].copy(),
            2 * (float(strength) * (1 - float(alpha)) * factors[retained].copy() + numerical_ridge))
        coefficients, objectives, criterion = solver.solve(
            design, response, datafit, penalty, w_init=coefficients.copy(), Xw_init=design @ coefficients)
        assert criterion <= solver.tol and np.isfinite(coefficients).all()
        full = np.zeros(original_X.shape[1])
        full[retained] = coefficients
        internal = full.copy()
        full[clinical] = np.linalg.solve(R, full[clinical] * np.sqrt(len(y))) - projection @ full[penalized]
        assert np.max(np.abs(original_X @ full - design @ coefficients)) < 1e-8
        residual = response - original_X @ full
        score = np.where(residual >= bandwidth, tau,
                np.where(residual >= 0, tau * residual / bandwidth,
                np.where(residual > -bandwidth, (1 - tau) * residual / bandwidth, tau - 1)))
        numerical_gradient = np.zeros_like(full)
        numerical_gradient[clinical] = numerical_ridge * R.T @ internal[clinical] / np.sqrt(len(y))
        numerical_gradient[penalized] = numerical_ridge * full[penalized] + projection.T @ numerical_gradient[clinical]
        gradient = -original_X.T @ score / len(y) + strength * (1 - alpha) * factors * full + numerical_gradient
        distance = np.maximum(np.abs(gradient) - strength * alpha * factors, 0)
        nonzero = np.abs(full) > 1e-9
        distance[nonzero] = np.abs(gradient[nonzero] + strength * alpha * factors[nonzero] * np.sign(full[nonzero]))
        maximum_kkt = float(distance.max())
        assert maximum_kkt < 1e-5
        value = datafit.value(response, full, original_X @ full) + strength * np.sum(
            factors * (alpha * np.abs(full) + (1 - alpha) * full**2 / 2)) + numerical_ridge * np.sum(internal**2) / 2
        beta.append(full)
        elapsed = time.perf_counter() - started
        diagnostics.append([index + 1, len(objectives), criterion, maximum_kkt, value, elapsed])
        print(f"smooth path {index + 1} KKT {maximum_kkt:.3g} seconds {elapsed:.3f}", flush=True)
    return {"beta": np.column_stack(beta), "diagnostics": np.asarray(diagnostics), "bandwidth": bandwidth}
