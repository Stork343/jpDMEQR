# QIF SEE

## Source

Indrabati Bhattacharya, Mohammad Alfrad Nobel Bhuiyan and Suneel Babu Chatla.
Automatic Variable Selection for Longitudinal Quantile Regression With Application
to Alzheimer's Disease Progression. Scandinavian Journal of Statistics 53,
1189–1205, 2026. DOI 10.1111/sjos.70077.

The complete 17-page publisher PDF was read through institutional access. Its
method is given by Equations 3–7, Algorithm 1 on page 1193 and the BIC on page
1194. The publisher lists no supplementary materials. The access record and
readable article text are retained in
`results/scratch/qif_see_fulltext_20261022/`.

## Statistical definition

For coordinate k and basis b, Equation 3 defines

`g_ikb = X_i[,k]' Gamma_i M_ib {Phi((y_i-X_i beta)/h) - (1-tau)}`.

`Gamma_i` is the diagonal conditional error density at zero. The article defines
this population quantity but supplies neither an estimator nor its numerical
assignment in the simulation section. The adapter therefore accepts explicit
observation-level `control$density_weights`. Formula tests use a genuine Normal
DGP with its known conditional density. Production fitting never substitutes a
truth-dependent weight or an invented density estimator.

The registered known basis consists of identity, off-diagonal ones and first
off-diagonal ones. The helper also exposes the source's independence and
exchangeable basis specifications. Each coordinate's B basis entries are
consecutive, as in Equation 3. The prescribed bandwidth is `n^(-1/2)`.

Algorithm 1 updates `delta_k=min(1,lambda/abs(beta_k)^(1+gamma))` at each iteration,
using `delta_k=1` when `beta_k=0`. It applies Equation 5 to the extended score,
`g_ikb^s=(1-delta_k)g_ikb+delta_k beta_k`. The sample mean and uncentred second
moment of this modified score define `Q=n G_s' C_s^(-1) G_s`. Within an iteration,
delta is held fixed in the approximate gradient `2 D' C_s^(-1) G_s` and approximate
Hessian `2 D' C_s^(-1) D`. The Newton update uses their ordinary inverses. Euclidean
step size determines convergence; iteration exhaustion is retained as failure.

The adapter obtains the initial coefficients from `rqPen::rq.pen.cv` with LASSO
and ten-fold cross-validation. Its returned intercept is recorded in the initial
object; the supplied slope design determines the coordinates subsequently fitted
by QIF. BIC is exactly `Q + df*log(n)`, with all declared lambda and gamma candidate
values retained. The source specifies a two-dimensional grid without numerical
endpoints. The project prospectively retains its eight lambda values from 0.05
through 1 and the declared gamma values, defaulting to 1. These project settings
are recorded separately from source formulas.

## Ordinary inverse domain

With B basis matrices and p coordinates, `C_s` is `Bp` by `Bp` and has rank at most
n. Every registered attempt with `Bp>n` records `qif_covariance_rank`, matrix
dimension, n, the rank upper bound and `rank_deficient_by_dimension`. Actual matrix
rank is recorded as unavailable when this algebraic check terminates fitting.
It applies for every possible conditional density input. For a computed modified
score matrix, numerical rank and the complete covariance are retained before any
ordinary inverse calculation. If a coordinate has `delta_k=1`, its B modified
score entries coincide; this can also make the covariance singular.

These domain failures remain attempted-method records. The source's ordinary
inverse and all registered method labels are preserved. A low-dimensional request
without density weights records `source_density_unspecified`.

## Verification

`test-benchmark-QIF-SEE.R` checks deterministic high-dimensional domain failure,
output fields, Equation 3 and Equation 5 ordering, finite-difference derivatives,
and the independence limit with negligible SEE weight. That limit is compared
with a separate convolution quantile optimizer on a genuine generated dataset.

`check_reference_fidelity_20261022.R` independently constructs all per-coordinate
basis scores and finite-difference derivatives, then computes Algorithm 1's
Newton step directly. It also verifies the algebraic rank boundary. The current
acceptance manifest is computed from these actual checks by
`accept_all_benchmarks_20261022.R`.

Allowed outputs are estimation, selection and fixed-effect prediction within the
source's computable domain. This adapter supplies no coordinate intervals.
