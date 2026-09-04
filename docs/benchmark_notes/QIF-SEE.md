# QIF-SEE — source note

## Reference

Bhattacharya, I., Bhuiyan, M.A.N., and Chatla, S.B. (2026). Automatic
Variable Selection for Longitudinal Quantile Regression With Application to
Alzheimer's Disease Progression. *Scandinavian Journal of Statistics*
53(3):1189-1205. DOI 10.1111/sjos.70077.

## Role in this project

Recent SJS longitudinal quantile variable-selection comparator. Estimation
and selection only; no coordinate coverage (the source paper supplies no
interval procedure for the selected coefficients).

## Implemented construction

The method combines the quadratic inference function (QIF) of Qu, Lindsay
and Li (2000) with the smooth-threshold estimating equations (SEE) of
Ueki (2009), as described in the source paper:

1. **QIF basis matrices** for the inverse working correlation (per cluster
   of size m): M1 = I_m (independence), M2 = off-diagonal ones
   (exchangeable), M3 = first off-diagonal ones (AR(1) lag-1).
2. **Induced-smoothed quantile score** per observation:
   `psi(r) = Phi(r/h) - (1 - tau)` with h = n^{-1/2} (Brown & Wang 2005
   smoothing; the source paper cites Brown & Wang and the companion QGEE
   literature uses the same convention). The analytic derivative
   `phi(r/h)/h` makes the estimating equations differentiable.
3. **Extended score** per cluster:
   `g_i(beta) = [ X_i' M1 psi_i ; X_i' M2 psi_i ; X_i' M3 psi_i ]`
   (3p x 1), stacked into G (3p x n).
4. **QIF objective** `Q_n(beta) = n gbar' Cbar^{-1} gbar` with
   `Cbar = n^{-1} sum_i g_i g_i'`.
5. **SEE weights** (Ueki 2009): `delta_k = min{1, lambda / |beta_tilde_k|^(1+gamma)}`
   with gamma = 1, applied to the QIF estimating equation U(beta):
   `(1 - delta_k) U_k(beta) + delta_k beta_k = 0`.
6. **Initial estimator** `beta_tilde` from `rqPen::rq.pen.cv`
   (lasso QR, 10-fold cross-validation), as stated in the source paper.
7. **Newton-Raphson iteration** alternating beta and Cbar with the analytic
   Jacobian `dU = n dgbar' Cbar^{-1} dgbar`.
8. **BIC tuning** over a lambda grid:
   `BIC = log(Q_n) + |M| log(n)/n` (|M| = number of nonzero coefficients).

## Formula-confirmation note (honest)

The Wiley full text was paywalled during implementation. The following
details are confirmed from indexed fragments and the cited predecessors:
QIF basis with independence/exchangeable/AR(1)/unstructured patterns;
SEE delta form `min{1, lambda/|beta_tilde_k|^(1+gamma)}`; modified score
`(1 - delta_j) g_ij + delta_j beta_j`; rqPen initial estimator with 10-fold
CV; Newton-Raphson-type update; BIC tuning selector. The following were
NOT confirmed from the full text and are implemented by standard choice:
(i) induced smoothing of the quantile score (adopted from Brown & Wang 2005,
   which the paper cites, and consistent with the companion QGEE method);
(ii) gamma = 1 (Ueki 2009 default); (iii) exact BIC constant (standard
   form used). If the full text becomes available, these three choices must
   be re-verified before final comparative claims.

## Allowed metrics

Estimation, selection, prediction. No coordinate coverage.

## Dependencies

- `rqPen` (initial estimator; CRAN);
- base R and project utilities for the QIF/SEE iteration.

## Acceptance evidence

`results/preflight/benchmarks/QIF-SEE/acceptance.json` from
`scripts/benchmarks/accept_qif_see.R` and unit test
`tests/testthat/test-benchmark-QIF-SEE.R`.

Fidelity checks:

1. Deterministic seed reproduces identical selected set and coefficients.
2. Active variables are retained and signs match the alternating design on
   a deterministic DGP.
3. Cluster-order permutation leaves results invariant up to numerical
   tolerance.
4. No coordinate coverage is produced (feasible = FALSE rows).
