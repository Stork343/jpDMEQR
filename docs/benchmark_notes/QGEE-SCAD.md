# QGEE-SCAD — source note

## Reference

Zu, T., Lian, H., Green, B., and Yu, Y. (2023). Ultra-High Dimensional
Quantile Regression for Longitudinal Data: An Application to Blood Pressure
Analysis. *Journal of the American Statistical Association* 118(541):97-108.
DOI 10.1080/01621459.2022.2128806.

## Implementation status

The authors maintain the reference R package **geeVerse** (CRAN 0.3.1,
GPL-3; GitHub zzz1990771/geeVerse). This adapter calls the official
implementation (`geeVerse::qpgee`) on the supplied training design, so the
fidelity status is `official_reference_implementation` rather than a
re-implementation. The package implements the paper's method:

- quantile penalised GEE with induced-smoothing score
  `Phi(sqrt(n) u / r_i) - (1 - tau)`, `r_i = ||X_i||` (row norm);
- SCAD penalty (a = 3.7) via LLA/MM with c = 1e-6 perturbation;
- working correlation structures: independence, exchangeable, AR(1), Tri,
  unstructured (moment-estimated from standardised quantile residuals);
- Newton-type updates with active-set shrinkage;
- HBIC or subject-level 5-fold CV lambda selection;
- SIS screening for ultra-high dimension (not invoked at the dimensions of
  this project's registry).

## Adapter decisions

- `corstr` defaults to `exchangeable` (the paper's primary setting); other
  structures are available through `control$corstr`.
- `method` defaults to `HBIC` (the paper's recommended selector).
- The design matrix carries no intercept column; the package derives its
  intercept flag from `attr(x, "assign")`, which the adapter sets explicitly
  to avoid the package's empty-attribute bug.
- Target interpretation: marginal longitudinal quantile estimating-equation
  target. Coverage against the profile target `beta_star` is permitted only
  in DGPs where the marginal and profile targets coincide, and the adapter
  itself never emits intervals (feasible = FALSE rows).
- A documented package-level detail: the working-correlation update uses
  `0.5 - 1{h < 1e-10}` for the residual sign in `geeVerse` internals, which
  differs from the paper's printed `tau - 1{...}` for tau != 0.5. The
  official package is treated as authoritative (fidelity to the maintained
  reference code).

## Allowed metrics

Estimation, selection, prediction; source-supported inference on its own
target only (not implemented in this adapter, so no coverage output).

## Dependencies

- `geeVerse` (CRAN), `quantreg` (transitively).

## Acceptance evidence

`results/preflight/benchmarks/QGEE-SCAD/acceptance.json` from
`scripts/benchmarks/accept_qgee_scad.R` and unit test
`tests/testthat/test-benchmark-QGEE-SCAD.R`.

Fidelity checks:

1. Deterministic seed reproduces identical output.
2. Active support recovered with correct alternating signs on deterministic
   data (marginal target; partial recovery accepted with recorded sign).
3. Working-correlation choice is recorded and both independence and
   exchangeable runs produce finite estimates.
4. No coordinate coverage is produced.
