# BIAS-ADJ-LQMM — source note

## Reference

Battagliola, M., Sørensen, H., Tolver, A., and Staicu, A.-M. (2022).
A bias adjusted estimator in quantile regression for clustered data.
*Econometrics and Statistics*, DOI 10.1016/j.ecosta.2021.07.003.
arXiv:2202.11501. Official R code: `bias_adj.zip` from
http://helle.sites.ku.dk/software/ (files `ran_int_bias_adj.R`,
`test_simulated_data.R`; archived in this repository's research material).

## Role in this project

Many-small-cluster direct comparator: bias-adjusted clustered quantile
regression. The source estimator conditions on estimated random effects
(centred BLP from LQMM) and corrects the resulting shrinkage bias with an
RW (resampling + wild) bootstrap. Used on a common low-dimensional screened
or oracle support only.

## Implemented algorithm (faithful to Algorithm 1 and official code)

1. **Step 1 (LQMM)**: fit a random-intercept linear quantile mixed model with
   `lqmm` (Gauss-Hermite quadrature, nK=15, type="normal",
   `lqmmControl(method="df", LP_max_iter=1000)`); extract the BLP random
   intercepts and centre them: `u_est = ranef - mean(ranef)`.
2. **Step 2 (rq)**: run standard quantile regression of
   `y - Z' u_est` on the fixed design at level `tau`; this is the two-step
   estimate `beta_hat` and its `nid` standard error `SE_obs`.
3. **Bootstrap** (B=100 by default, as in the paper):
   - wild weights `w ~ {2(1-tau) w.p. 1-tau; -2tau w.p. tau}` (zero
     tau-quantile);
   - cluster random effects resampled with replacement from the centred BLP
     distribution `{u_est_i}`;
   - bootstrap response `Y* = X'beta_hat + Z' u* + w|residual|`;
   - *oracle replicate*: rq of `Y* - Z'u*` on X (no LQMM);
   - *two-step replicate*: re-run the full LQMM+rq pipeline on `Y*`.
4. **Bias adjustment**: `beta_adj = 2 beta_hat - mean(beta_two_step*)`
   (Efron-Tibshirani bootstrap bias correction).
5. **SE adjustment**: `SE_adj = SE_obs * SD(two-step*) / SD(oracle*)`.
6. **Intervals**: SE-adjusted Wald interval `beta_adj +/- z SE_adj`
   (recommended in the paper) and basic bootstrap interval
   `(2 beta_hat - q_{1-a/2}, 2 beta_hat - q_{a/2})`.

## Deviations / decisions

- Random-intercept-only support (q=1 nuisance column), matching the source
  `ran_int_bias_adj.R`. A random-slope extension is explicitly refused with
  a `failed` status rather than silently approximated.
- The fixed-effect formula is `y ~ 0 + X` (no forced intercept) so the
  coefficient names align with the supplied design; the source code's
  intercept handling is reproduced when the design contains one.
- Bootstrap count is B=100 unless overridden by `control$B`; a run with fewer
  than 50% successful bootstrap replications returns `failed` (faithful
  failure reporting rather than a degraded estimate).
- Tuning uses training clusters only; seeds are deterministic.

## Allowed metrics

Low-dimensional estimation, bias, prediction where defined, and coverage
only through the faithful SE-adjusted/basic-bootstrap interval procedure.

## Dependencies

- `lqmm` (random-intercept LQMM step);
- `quantreg` (second-step and bootstrap rq fits).

## Acceptance evidence

`results/preflight/benchmarks/BIAS-ADJ-LQMM/acceptance.json` from
`scripts/benchmarks/accept_bias_adj_lqmm.R` and unit test
`tests/testthat/test-benchmark-BIAS-ADJ-LQMM.R`.

Fidelity checks:

1. Two-step estimate matches a manual LQMM+rq pipeline on a fixed seed.
2. With negligible random-effect variance, the bias adjustment has a small
   effect (limits to ordinary QR with a small cluster correction).
3. Zero-quantile property of the wild weights is verified empirically.
4. Deterministic seed reproduces identical output.
5. Failure is returned (not a degraded estimate) when bootstrap success
   fraction is too low.
