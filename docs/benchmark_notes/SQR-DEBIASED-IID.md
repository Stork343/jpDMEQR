# Smoothed quantile regression inference

`SQR-DEBIASED-IID` implements the independent-observation procedure of
Yan, Wang and Zhang (2023), *Confidence Intervals and Hypothesis Testing
for High-dimensional Quantile Regression: Convolution Smoothing and
Debiasing*, Journal of Machine Learning Research 24(245), 1–49.
The [journal record](https://www.jmlr.org/papers/v24/22-1217.html) supplies
the published article.

## Statistical definition

The Gaussian-smoothed check loss is
`u * (tau - Phi(-u/h)) + h * phi(u/h)`. Its L1-penalised empirical
minimiser uses an unpenalised working intercept. With N observations,
the score and Hessian are

\[
 g=N^{-1}\sum_i[\Phi\{(X_i^T\widehat\beta-Y_i)/h\}-\tau]X_i,
 \qquad
 H=(Nh)^{-1}\sum_i\phi\{(X_i^T\widehat\beta-Y_i)/h\}X_iX_i^T.
\]

For a common tolerance gamma, every required CLIME direction minimises
its L1 norm subject to `max(abs(H %*% omega - e_k)) <= gamma`.
Pairwise minimum-magnitude symmetrisation gives the matrix W used in
the coefficient correction. Its requested row w gives

\[
 \widetilde\beta_k=\widehat\beta_k-w^Tg,\qquad
 \widehat{\operatorname{se}}_k^2
 =\tau(1-\tau)w^T\widehat\Sigma w/N,
 \qquad \widehat\Sigma=N^{-1}\sum_iX_iX_i^T.
\]

The bandwidth is
`max(sqrt(tau*(1-tau))*sqrt(log(p))*N^(-3/10), 0.05)`, with p the
working design dimension. The initial penalty is
`c_lambda*sqrt(tau*(1-tau)*log(p)/N)`.

## Project tuning protocol

Five whole-subject folds choose the initial penalty multiplier using
equal-subject pinball loss. Five subject folds choose the common
precision tolerance by mean squared full-matrix inverse error,
`sum((H_validation %*% W - I)^2)/p`. Each precision-training fold
reselects its initial penalty within its own subjects.

The multiplier grids and fold-specific scales come from the frozen
registry. The candidate records retain the configured values and
their actual tolerances. Values at or above one have the recorded
zero-direction exclusion. The selected tolerance is then used for
the complete full-training precision matrix and coordinate intervals.

## Numerical implementation

Matrix-free FISTA computes the initial estimator. HiGHS solves the
CLIME programmes through a compiled active-row and active-column
driver, with full primal-dual verification. Certified column-loss
bounds permit early elimination of candidates that cannot minimise
the complete validation criterion. Eligible candidates and the final
fit retain their complete symmetrised matrices.

The bound, its floating-point allowance and the retained first-direction
path budget are documented in
[Certified computation for SQR](../SQR_COMPUTATION_20261030.md).
The execution specification is
`config/submission_comparator_numerics_20261030.yml`.

## Outputs and validation

The estimator reports estimation and selection metrics, prediction,
and intervals with the `naive_iid` variance identity. Each scenario
retains its declared target. Predictor-dependent scale settings use
the separately calculated marginal quantile projection.

Validation compares complete saved datasets, including the chosen
penalty and tolerance, coefficients, standard errors and endpoints.
Independent R HiGHS calculations check the compiled CLIME solutions.
OSQP checks the validation bound against box-constrained quadratic
programmes. A registered 2000-predictor case independently verifies
all selected-fold losses and the final coefficient correction.

Required libraries are RcppEigen, the R HiGHS interface and the compiled
HiGHS 1.15.1 library. The retained parametric path uses fastclime
1.4.1.1. `scripts/00_install_dependencies.R` declares the installation
steps, and numerical checks are saved under
`results/preflight/sqr_acceleration_20261030/`.
