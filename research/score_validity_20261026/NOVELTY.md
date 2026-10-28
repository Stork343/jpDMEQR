# Method Contributions and Related Work

The contribution is a computational and inferential construction for the complete quantile profile with many independent small clusters. The working nuisance ridge changes the profiled score geometry, and the selected model must approximate the full population target. These two features determine the algorithm and its theory.

## Published foundations

| Construction | Published foundation | Role in this method |
| --- | --- | --- |
| Orthogonal variable search | Ing and Lai develop OGA and model selection for high dimensional sparse linear regression. Cha, Chiang and Sasaki combine OGA with high dimensional model selection and debiased inference. | The selected subproblem is a complete clustered quantile profile. Its dual score drives the search, and a prediction ridge supplies an extension bound and an optimisation certificate. |
| Outcome and adjustment selection | Belloni, Chernozhukov and Hansen develop post double selection with inference under imperfect control selection. | The adjustment set combines a profile curvature graph with conditional covariate regression. The latter can preserve a weak coefficient's target even when prediction omits other coefficients. |
| Smoothed quantile inference | Yan, Wang and Zhang construct debiased smoothed quantile estimates and confidence procedures in high dimension. | Repeated observation profiling can produce a term linear in bandwidth. Two complete convex profile fits remove that term, with their joint cluster covariance retained. |
| Primal dual optimisation | Malitsky and Pock provide the linesearch algorithm and its convergence analysis. | The complete joint quantile objective has explicit proximal maps and a computable primal dual gap. The statistical algorithm uses these existing optimisation results. |

Ing and Lai's paper is *A Stepwise Regression Method and Consistent Model Selection for High-Dimensional Sparse Linear Models*, Statistica Sinica 21, 1473–1513, DOI 10.5705/ss.2010.081. The author institution's [publication record](https://scholars.nthu.edu.tw/esploro/outputs/journalArticle/A-stepwise-regression-method-and-consistent/9957776627806774) gives its scope and bibliographic details. The combined OGA inference construction is described by [Cha, Chiang and Sasaki](https://arxiv.org/abs/2108.09520).

[Belloni, Chernozhukov and Hansen](https://arxiv.org/abs/1201.0224) establish the use of outcome and adjustment selection for inference when control selection is imperfect. Their publication is The Review of Economic Studies 81, 608–650, DOI 10.1093/restud/rdt044. The quantile smoothing comparison is [Yan, Wang and Zhang](https://www.jmlr.org/papers/v24/22-1217.html), Journal of Machine Learning Research 24(245), 1–49, 2023. The optimisation source is [Malitsky and Pock](https://epubs.siam.org/doi/abs/10.1137/16M1092015), SIAM Journal on Optimization 28, 411–432, 2018.

## Complete profile computation

The cluster nuisance minimisation has a finite solution construction for an intercept and for an intercept with a time slope. Its dual certificate is checked at each accepted solution. The joint selected quantile problem supplies an additional certificate without repeatedly solving a complete nuisance profile inside every high dimensional optimisation step. These constructions are derived in `positive_curvature_20261025/FINITE_PROFILE.md` and `ORTHOGONAL_METHOD.md`.

The inferential calculation uses small selected positive Hessians. Both bandwidth problems remain convex throughout optimisation. Each target's variance then combines the two fitted cluster contributions, including their covariance. Computational cost can be traced to full design scans during selection and selected matrix operations during inference; `RESULTS.md` reports the actual local timings.

## Preserving the full target

Prediction accuracy alone can omit a covariate needed to adjust a target coefficient. The method also uses conditional covariate associations, while retaining the profile graph approximation for the general risk geometry. The proof controls the selected model target by the product of outcome and contrast approximation errors, plus a quadratic outcome error.

For independent Gaussian covariate rows, the selected population profile coefficient admits an explicit conditional projection formula. Inclusion of the target's conditional neighbours makes its target coordinate equal to the full coefficient at every bandwidth. `JOINT_THEORY.md` proves that statement and gives the primitive learning bounds leading to asymptotic normality. The result does not require a minimum nonzero outcome coefficient.

## Contribution wording for the manuscript

We develop feasible coordinate inference for high dimensional quantile profiles with many small clusters. Conditional covariate selection controls the displacement of the population target after dimension reduction, and paired positive curvature fits remove the first order profile smoothing term. The resulting estimator has a cluster influence representation under explicit sparse learning and density conditions. A finite nuisance solver and certified joint optimisation make the method practical on a local computer.

The claims are attached to the complete profile construction, target preservation and the derived statistical guarantees. OGA, double selection and primal dual linesearch retain the attributions above.
