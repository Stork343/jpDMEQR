# Yeast quantile prediction comparisons

The data are observed replicate growth measurements from Bloom et al. (2013), with one segregant per tetrad. All prediction models use the same frozen strain roles. A strain receives weight one, and its observations share that weight equally. Training transformations and tuning use whole-strain folds. Calibration estimates replicate-specific residual quantiles from out-of-fold training predictions. Evaluation outcomes are used only to calculate prediction losses.

## Exact quantile Lasso

For N rows and n strains, the fitted criterion is

\[
 n^{-1}\sum_i m_i^{-1}\sum_j\rho_\tau(y_{ij}-x_{ij}^{T}\beta)
 +\lambda\sum_k a_k|\beta_k|.
\]

The nonnegative homogeneity of check loss permits the existing native HiGHS solver to use rows and responses multiplied by N/(n m_i). This preserves the displayed criterion, including its penalty. Each internal training fold estimates its own maximum penalty. Forty common ratios from 1 to 0.01 define the search path. Three whole-strain folds select the ratio by equal-strain pinball loss. The unpenalized adjustment columns and complete genomic and expression blocks follow the corresponding application configuration.

All-coordinate dual feasibility and the original primal-dual gap are checked at every path point. `check_pooled.R` compares the weighted objective with `quantreg::rq.fit.lasso` on observed data and checks the zero-penalty limit with `quantreg::rq.fit.fnb`. The numerical engine and its earlier complete-path checks are described in `SOYNAM-POOLED-LASSO.md`. This application adapter provides prediction comparisons; it does not construct confidence intervals.

## Quantile forests

The forest comparator directly calls `grf::quantile_forest` version 2.6.1. Its source documentation is https://grf-labs.github.io/grf/reference/quantile_forest.html and the statistical reference is Athey, Tibshirani and Wager (2019), *Generalized random forests*, Annals of Statistics 47, 1148–1178, DOI 10.1214/18-AOS1709.

The native `clusters` and `equalize.cluster.weights=TRUE` options handle repeated observations. The grid combines quantile and regression splitting, mtry values proportional to sqrt(p)+20, minimum leaf sizes 5 and 10, and honest or non-honest trees. Each fit has 2000 trees. Three whole-strain training folds select settings for each of the three quantiles. All candidate settings and final choices are retained. Forest predictions use the same training residual calibration and equal-strain evaluation as the linear methods.

`check_forest.R` verifies deterministic predictions and ordered native quantiles on observed training strains. `report_comparators.R` reconstructs predictions from saved forests, verifies their row identities and losses, and constructs paired strain-bootstrap comparisons. The separately calibrated quantiles are evaluated individually. No forest interval-coverage claim is made.

## Interpretation of the comparisons

The tables are exploratory comparisons from the frozen discovery roles. Selection among culture conditions is recorded in `yeast_screening` and `quantile_comparisons_20261101`. Their pointwise bootstrap intervals do not certify a condition selected from the collection. Formal validation uses the separately retained strains and the registered sequential testing allocation after the scientific question and complete candidate order have been fixed.
