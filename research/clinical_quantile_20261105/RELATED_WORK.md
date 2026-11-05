# Closest method comparison

## Longitudinal penalised quantile estimation

Zu, Lian, Green and Yu, Ultra-High Dimensional Quantile Regression for Longitudinal Data An Application to Blood Pressure Analysis, JASA 118, 97–108, 2023, DOI 10.1080/01621459.2022.2128806. The repository has an audited official geeVerse adapter and a method note at `docs/benchmark_notes/QGEE-SCAD.md`. This work already combines high dimensional variable selection with quantile estimating equations and working correlation. Correlation weighting alone therefore cannot identify the contribution of the new construction.

Song, Han, Fu and Wang, Penalized weighted smoothed quantile regression for high dimensional longitudinal data, Statistics in Medicine 43, 2007–2042, 2024, DOI 10.1002/sim.10056. The publisher abstract describes penalised convolution smoothing and a two step weighted estimator that incorporates within patient dependence. Formula level comparison requires its full estimation and inference construction; the abstract establishes that weighted smoothed longitudinal quantile estimation is already studied.

## Orthogonal quantile estimation

Qiu, Zhu, Wang and Tang, Heterogeneous Quantile Treatment Effect Estimation for Longitudinal Data with High Dimensional Confounding, arXiv:2508.16326v2, 31 August 2026. Section 2.1 uses a location shift condition and states that only the intercept varies across quantiles. Equations 3 and 4 use a treatment residual in the quantile score and average within each patient. Forest weights localise the treatment effect by baseline modifiers. Section 2.4 uses Bootstrap of Little Bags for uncertainty.

The current construction studies a fixed coefficient at each quantile, permits varying conditional densities, and residualises its target through F_i C_i^(-1) F_i. It includes the full patient covariance in the instrument and uses patient influence contributions for uncertainty. Establishing a distinct contribution requires a primitive rate theorem for those estimated density and covariance terms, alongside comparisons with the longitudinal estimating equation methods.

Giessing and Wang, Debiased Inference on Heterogeneous Quantile Treatment Effects with Regression Rank Scores, JRSSB 85, 1561–1588, 2023, arXiv:2102.01753v5. Its abstract establishes an L1 regression adjustment, rank score correction, weak convergence and semiparametric efficiency. A new efficiency claim must be compared at the level of the target and instrument class, including the corresponding restrictions on the conditional distribution.

## Computation

Lavancier and Rochet, A general procedure to combine estimators, Computational Statistics and Data Analysis 94, 175–192, 2016, arXiv:1401.6371. The authors use an estimated mean squared error matrix to combine estimators of a common parameter and establish asymptotic optimality. The convex aggregation in the present construction uses the patient influence covariance after the two component estimators have a centred joint expansion. The aggregation principle is inherited from this literature.

Yi and Huang, Semismooth Newton Coordinate Descent Algorithm for Elastic Net Penalized Huber Loss Regression and Quantile Regression, JCGS, DOI 10.1080/10618600.2016.1256816. Penalised quantile paths are established numerical objects. The current certified computation is an implementation contribution only if its full objective, accuracy and end to end costs improve on relevant implementations at the required dimensions. Timing one favourable path does not establish that comparison.

The reviewed sources identify existing components and the comparisons still required. They do not establish a priority claim for the current combination.
