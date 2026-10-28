# PROFILE-DQR

The procedure is defined in `docs/METHOD_SPECIFICATION.md` and the manuscript
sections on profile estimation and coordinate inference. Its population
parameter minimises the unsmoothed complete ridge-profile risk. The working
matrix Lambda remains fixed.

`fit_benchmark_profile_pair_v2` calls the public convex selection and refit
API. The initial penalty is the sum of the cluster stochastic contribution
and `mean_score_multiplier * h_est`, multiplied by the fixed final score
loadings and the prospective penalty multiplier. The two loading passes
and both penalty contributions remain in the output. The inference score
is `2 * g_half - g_full`, with separate nuisance optima at one coefficient
vector. Its Jacobian uses the same combination of exact Schur Hessians.

Precision directions solve the row-wise L1 programme over the requested
frozen multiplier grid. Candidates outside `(0,1)` are recorded before
cluster-fold validation. The update and centred cluster sandwich use the
same combined score. CLARABEL numerical controls and the achieved inverse
defects are retained.

Allowed metrics are estimation, selection, prediction, prespecified-coordinate
coverage, theory diagnostics and computation. The source code is part of
this project. R, Matrix and clarabel provide the numerical operations;
their package licenses apply to those dependencies.

`tests/testthat/test-benchmark-profile-pair.R` checks deterministic public-API
mapping, full-design predictions including permuted columns, independent
profile recomputation and the same-score correction and variance identities.
Its high-dimensional example uses a genuine random intercept and time slope
with a non-diagonal Lambda. The paired geometry and population tests supply
the differential and limiting calculations. The acceptance script is
`scripts/benchmarks/accept_profile_pair_20261022.R`.
