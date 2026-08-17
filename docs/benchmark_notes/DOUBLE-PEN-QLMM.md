# DOUBLE-PEN-QLMM — source note

## Reference

Li, Y., Liu, Y., and Luo, Y. (2020). Variable selection for joint fixed and
random effects in quantile mixed models. *Journal of Systems Science and
Complexity* 33:2080-2102. DOI 10.1007/s11424-020-9065-4.

## Role in this project

Direct penalised mixed-quantile comparator: jointly selects fixed and random
effects with two L1 penalties. Estimation/selection comparator; no coverage
claim by default (the source paper supplies no interval procedure).

## Implemented algorithm (faithful to paper eq. (12) and Section 3.1)

Objective:
```
L(beta, alpha) = sum_i sum_j rho_tau(y_ij - x_ij'beta - z_ij'alpha_i)
                 + lambda_beta sum_l |beta_l|
                 + lambda_alpha sum_i sum_t |alpha_it|
```
Two-block alternating L1 quantile regression:

1. Initialise `alpha = 0`; fit beta^(0) by L1 quantile regression of y on X.
2. Iterate until `mean_k |beta^(m+1) - beta^(m)| < eps` (paper uses 1e-4):
   - **alpha-step**: for each cluster, L1 quantile regression of
     `r_ij = y_ij - x_ij'beta` on the cluster's z_ij (cluster-separable);
   - **beta-step**: L1 quantile regression of `y_ij - z_ij'alpha_i` on X.
3. Coefficients below 1e-6 are set exactly to zero (paper threshold).

Tuning (paper eq. (13), one lambda at a time):

- SIC(lambda) = log(S_M/N) + log(N)/(2N) * |M|
- GACV(lambda) = S_M / (N - |M|)
- `S_M` = sum of check losses at the fit; `|M|` = number of nonzero
  coefficients in the current block; grid = deterministic positive values.
- Implementation selects `lambda_beta` by SIC on a frozen grid with
  `lambda_alpha` fixed at the grid midpoint, mirroring the paper's
  one-at-a-time rule; the selected values are recorded.

## Deviations / decisions

- No official code exists for the paper (verified: no repository, no
  supplementary code). Implementation is written from the printed formulas.
- The paper's per-block L1 quantile regression is solved with
  `quantreg::rq.fit.lasso`; zero-penalty blocks fall back to plain
  `quantreg::rq.fit`.
- The paper's simulation threshold `|coef| < 1e-6 => 0` is reproduced.
- Stopping monitors beta only (as in the paper); alpha is updated inside the
  same iteration.
- Random effects are estimated as penalised coefficients (paper's
  "random-effect estimation = penalised least squares" view); no variance
  component D is estimated, matching the paper.
- No inferential output: `beta_tilde`, `se`, and intervals are NULL/NA, and
  `inference_object$table$feasible` is FALSE for all rows, so the schema
  runner records no coverage rows.

## Allowed metrics

Estimation, selection, prediction. No coordinate coverage.

## Dependencies

- `quantreg` (`rq.fit.lasso`, `rq.fit`).

## Acceptance evidence

`results/preflight/benchmarks/DOUBLE-PEN-QLMM/acceptance.json` from
`scripts/benchmarks/accept_double_pen_qlmm.R` and unit test
`tests/testthat/test-benchmark-DOUBLE-PEN-QLMM.R`.

Fidelity checks:

1. Objective matches paper eq. (12) exactly on deterministic data.
2. Zero-penalty limit behaves as the unpenalised alternating estimator.
3. Active-coefficient signs recovered on deterministic DGP.
4. Deterministic seed reproduces identical output.
5. No inferential claim is produced (feasible = FALSE rows).
