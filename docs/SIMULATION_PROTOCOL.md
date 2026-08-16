# Preregistered simulation protocol

## 1. Why the simulation study is necessary

The proposed method combines high-dimensional regularisation, low-dimensional cluster-specific nuisance effects, convolution smoothing, profile curvature estimation, sparse inverse-information estimation and cluster-level uncertainty quantification. Asymptotic normality alone does not show that these components work jointly at sample sizes encountered in longitudinal studies. The simulation study therefore has five distinct scientific functions.

First, it must demonstrate the **necessity of profiling**. A pooled quantile Lasso can be robust to marginal outliers yet still be biased or inefficient when persistent subject heterogeneity is ignored. The study must identify when the profile adjustment materially changes estimation, prediction and inference.

Second, it must provide a **finite-sample audit of the proof architecture**. The paper relies on target approximation, penalised estimation, effective-Hessian concentration, sparse precision recovery, a Bahadur representation and a cluster sandwich variance. Each step has a numerical diagnostic. A table of prediction errors alone cannot validate the theoretical mechanism.

Third, it must compare the method with **classical, direct recent and oracle procedures**. The relevant question is not whether the method beats a deliberately weak pooled baseline in every setting. The study must show where it gains from profile geometry, where a marginal longitudinal method is competitive, and how far it lies from oracle performance.

Fourth, it must evaluate the **error-distribution robustness expected of quantile regression**. Gaussian, heavy-tailed, sharply peaked, skewed, contaminated and heteroskedastic errors are required. Every distribution is shifted so its conditional target quantile is exactly zero.

Fifth, it must expose **failure boundaries**. Dense inverse information, strong collinearity, nuisance misspecification, random slopes, weak signals and extreme quantiles may challenge the method. Failures and numerical instability are reported rather than hidden.

## 2. General principles

- Independent unit: cluster/subject.
- Primary sample-size symbol: `n`, the number of clusters.
- Cluster sizes are bounded in all primary experiments.
- Common random numbers are used for paired method comparisons.
- The true active set is known only to the data generator and oracle methods.
- All non-oracle tuning uses training clusters only.
- No method is tuned by empirical coverage or empirical bias.
- Failed fits remain in the denominator and are assigned a failure class.
- The raw replication-level output is immutable; all summaries are regenerated from it.

## 3. Baseline data-generating mechanism

For cluster `i=1,...,n` and observation `j=1,...,m_i`, generate

\[
Y_{ij}=X_{ij}^{\mathsf T}\beta^0+Z_{ij}^{\mathsf T}b_i+\varepsilon_{ij}.
\]

### 3.1 Cluster sizes

The baseline uses

```text
m_i ~ discrete Uniform{3,4,5,6,7,8}
```

independently across clusters. Sensitivity settings use `{2,3,4}` and `{6,...,12}` while remaining bounded. An additional imbalance setting draws from a truncated shifted geometric distribution and records the realized distribution.

### 3.2 Fixed-effect design

Unless a module states otherwise,

\[
X_{ij}\sim N_p(0,\Sigma_X),\qquad
(\Sigma_X)_{kl}=\rho_X^{|k-l|},\quad \rho_X=0.5.
\]

The support is `S0={1,...,s}`. Nonzero coefficients alternate in sign:

\[
\beta_k^0=a_s(-1)^{k+1},\qquad k\in S_0.
\]

The default signal is `a_s=0.75`; weak and strong settings use `0.40` and `1.10`. Predictors are centred and scaled within the data-generating population. No penalised fixed intercept is included, which helps ensure that nuisance shrinkage does not alter slope targets in the correctly specified baseline.

Alternative designs:

- independent: `rho_X=0`;
- strong AR(1): `rho_X=0.8`;
- block correlation: blocks of 25 variables with within-block correlation 0.7;
- approximate-factor design: three latent factors plus idiosyncratic noise;
- sparse-precision design: Gaussian graphical-model covariance with banded inverse;
- dense-precision stress design: equicorrelated or factor covariance producing non-sparse inverse rows.

The dense-precision design is a stress test, not a setting in which the Dantzig assumptions are expected to hold.

### 3.3 Random-effect design

Random-intercept setting:

\[
Z_{ij}=1,\qquad b_i\sim N(0,\sigma_b^2).
\]

Random-intercept--slope setting:

\[
Z_{ij}=(1,t_{ij})^{\mathsf T},
\]

where visit time is generated independently of `X`, centred within the population and scaled to unit variance. The random effects have covariance

\[
\Psi_b=
\begin{pmatrix}
\sigma_0^2 & \rho_b\sigma_0\sigma_1\\
\rho_b\sigma_0\sigma_1 & \sigma_1^2
\end{pmatrix}.
\]

Defaults are `sigma0=1`, `sigma1=0.5` and `rho_b=0.4`. Variance-strength settings multiply `Psi_b` by `0.25`, `1` and `4`.

Non-Gaussian random effects include:

- multivariate `t_5` rescaled to covariance `Psi_b`;
- symmetric two-component normal mixture with the same covariance;
- skewed random intercept obtained from a centred log-normal distribution.

### 3.4 Error distributions

For every target quantile `tau`, an error generator returns `epsilon = U - F_U^{-1}(tau)`. When heteroskedastic scaling is applied, the scale is strictly positive, preserving conditional `tau`-quantile zero.

Required distributions:

1. **Gaussian**: standard normal.
2. **Student t3**: rescaled to unit variance.
3. **Laplace**: double exponential rescaled to unit variance.
4. **Skewed chi-square**: standardized `chi-square_3`.
5. **Contaminated normal**: `0.9 N(0,1)+0.1 N(0,25)`, centred at the numerical mixture quantile.
6. **Asymmetric Laplace**: generated with zero `tau`-quantile and rescaled to unit standard deviation.
7. **Heteroskedastic Gaussian**: `epsilon_ij=exp(0.35 X_{ij,1}) U_ij`, with `U` quantile-centred.
8. **Heteroskedastic t3**: same positive scale with t3 innovations.
9. **Within-cluster copula dependence**: Gaussian copula with AR(1) parameter 0.4 and a selected marginal distribution, used only in a dependence-misspecification sensitivity experiment.

## 4. Target definition in simulations

### 4.1 Correctly specified slope scenarios

For baseline designs with centred `X` independent of `Z`, `b_i` and the quantile-centred error, evaluation uses `beta0` for fixed-effect slopes. The intercept is not a target coordinate.

### 4.2 Misspecified scenarios

The following scenarios target `beta_star`, not automatically `beta0`:

- omitted random slope;
- random scale driven by the random intercept;
- correlation between `X` and random effects;
- informative cluster size;
- deliberately misspecified nuisance penalty/direction;
- nonlinear random-effect contribution.

For these settings, `beta_star` is approximated by a population Monte Carlo optimisation using at least 100,000 independently generated clusters and a low-dimensional target submodel containing all active variables plus designated null controls. The population run is repeated with independent seeds; the difference between repeats must be negligible relative to the reported Monte Carlo standard error. The approximation file and optimisation diagnostics are saved.

## 5. Methods

### 5.1 Proposed method

`PROFILE-DQR`: corrected ridge-profile estimator, row-wise Dantzig direction, one-step correction and cluster sandwich variance as specified in `docs/METHOD_SPECIFICATION.md`.

Variants used for diagnosis:

- `PROFILE-DQR-TRUE-SUPPORT`: proposed method restricted to the true support plus target null coordinates;
- `PROFILE-DQR-TRUE-NUISANCE`: uses generated random effects in the score/curvature to isolate nuisance-estimation cost;
- `PROFILE-DQR-POP-H`: uses a high-accuracy population effective Hessian/direction;
- `PROFILE-DQR-UNSMOOTHED-SCORE`: optional sensitivity using the nonsmooth score in the final cluster variance while retaining smoothed curvature;
- `PROFILE-DQR-SPLIT`: subject-split screening followed by inference on independent clusters.

### 5.2 Classical methods

1. `POOLED-QR-LASSO`: l1-penalised quantile regression ignoring clusters.
2. `LQMM`: linear quantile mixed model using the `lqmm` package, fitted on the true support or a common screened set because it is not designed for `p >> n`.
3. `BIAS-ADJ-LQMM`: many-small-cluster bias-adjusted two-step/bootstrap procedure where computationally feasible.
4. `DOUBLE-PEN-QLMM`: double-penalised quantile mixed-effects estimator of Li, Liu and Luo (2020), reimplemented from the paper if no maintained reference code is available.

### 5.3 Recent high-dimensional/longitudinal methods

1. `QGEE-SCAD`: ultra-high-dimensional longitudinal quantile penalised GEE/SCAD method of Zu, Lian, Green and Yu (2023), with its working-correlation choices.
2. `SQR-DEBIASED-IID`: debiased convolution-smoothed high-dimensional quantile regression of Yan, Wang and Zhang (2023), applied naively at the observation level; this isolates the cost of ignoring clustering while using a modern inferential method.
3. `QIF-SEE`: automatic longitudinal quantile variable-selection method of Bhattacharya, Bhuiyan and Chatla (2026), included for estimation/selection when an implementation is reproducible. It is not treated as a coordinate-wise debiased-inference comparator unless its inferential target matches.

### 5.4 Predictive supplemental method

`BAYES-MIXED-LASSO` may be included as a predictive shrinkage comparator. Because it is generally a Gaussian mean model rather than a quantile-inference method, it is not a primary coverage comparator.

### 5.5 Fairness rules

- Methods receive the same generated data and, where applicable, the same training/validation subject split.
- Methods unable to fit the full dimension receive the same prespecified screened set; this limitation is labelled.
- Oracle methods are visibly labelled and excluded from claims of practical superiority.
- A method is not removed because it performs well or poorly.
- If an implementation is unavailable, the benchmark is marked `implementation_required`; it is not silently replaced by a different method with a similar name.

## 6. Simulation modules

The final registry is modular rather than a prohibitively large full factorial design.

### Module A: finite-sample theory scaling

Purpose: verify convergence trends and the link between the theorem and finite samples.

```text
n in {100, 200, 400}
p in {500, 2000}
s in {5, 10}
tau = 0.50
q = 1
error = Gaussian
rho_X = 0.50
m_i in {3,...,8}
```

Core methods: proposed, pooled QR-Lasso, true-support proposed, population-H proposed, SQR-debiased-iid.

Diagnostics: `l1/l2` error, Hessian max error, Dantzig residual, precision-row error, Bahadur remainder, empirical SD, mean SE, coverage, interval length and log-log rate slopes.

Core configurations use `B=1000`. High-cost `p=2000` precision diagnostics may use `B=500` if documented before the run.

### Module B: quantile and error robustness

Purpose: establish the usual robustness profile expected of quantile methods.

```text
n = 200
p = 1000
s = 8
tau in {0.25, 0.50, 0.75}
error in {Gaussian, t3, Laplace, skewed-chi2, contaminated-normal,
          asymmetric-Laplace, heteroskedastic-Gaussian, heteroskedastic-t3}
q = 1
```

Use a balanced fractional registry rather than all 24 combinations when runtime is limiting, but every error family and every target quantile must appear at least twice. Core error settings use `B=500`; Gaussian/t3/contaminated settings at all three quantiles use `B=1000`.

### Module C: nuisance structure and many-small-cluster difficulty

Purpose: test the central scientific problem.

Settings:

1. random intercept, low/medium/high variance;
2. correlated random intercept and time slope;
3. t5 random effects;
4. mixture random effects;
5. omitted random slope in the fitted model;
6. working `Lambda` misspecification;
7. small clusters `{2,3,4}`;
8. moderate bounded clusters `{6,...,12}`.

Primary comparison includes LQMM, bias-adjusted LQMM, double-penalised QLMM, QGEE-SCAD, proposed and oracle proposed.

### Module D: design and precision stress

Purpose: test sparse inverse-information assumptions.

```text
rho_X in {0, 0.5, 0.8}
design in {AR1, block, sparse-precision, factor, dense-precision}
```

The dense-precision setting is expected to challenge Dantzig estimation. Report constraint feasibility, row norms, interval length and coverage degradation. Do not describe failure under assumption violation as an implementation bug.

### Module E: smoothing and tuning sensitivity

Purpose: separate bandwidth bias, curvature variance and numerical conditioning.

```text
c_h in {0.75, 1.00, 1.25}
Lambda multiplier in {0.25, 0.5, 1, 2, 4}
c_lambda in {0.25, 0.5, 1, 2}
c_mu in {0.5, 1, 2, 4}
```

Only one dimension is varied at a time around the baseline. Report coverage with Monte Carlo uncertainty, but tuning choices are never made by choosing the value closest to nominal coverage.

### Module F: screening and ultra-high dimension

Purpose: evaluate computational dimension reduction without making screening the main contribution.

```text
n in {200, 400}
p in {5000, 10000}
s = 10
tau in {0.25,0.5,0.75}
```

Compare:

- no screening where computationally feasible;
- outcome-blind variance filter;
- cluster-adjusted IQR-SIS as a heuristic;
- independent subject-split quantile screening;
- oracle support.

Report sure-inclusion probability, active-set loss, selected dimension, prediction and inference on held-out subjects. Ordinary same-sample Wald intervals after heuristic screening are not labelled confirmatory.

### Module G: computational scaling

Purpose: provide realistic runtime and memory information.

Vary `n`, `p`, `q`, selected dimension, number of target coordinates and number of worker processes. Use `B=50` or repeated timed runs. Report median, IQR and maximum runtime, peak memory, iterations and failure rates. Do not combine these runs with statistical performance tables.

## 7. Replication counts and Monte Carlo uncertainty

- geometry tests: deterministic grid over at least 20 seeds;
- pilot registry: `B=200` per configuration;
- core inference scenarios: `B=1000`;
- secondary robustness scenarios: `B=500`;
- computational scaling: at least 20 timed repetitions;
- expensive population-target approximation: independently repeated twice.

For a coverage estimate `p_hat`, report

\[
MCSE=\sqrt{p_hat(1-p_hat)/B}
\]

and a binomial confidence interval. For means, report the Monte Carlo standard error `sd/sqrt(B)`. Paired method differences use the replication-level paired standard deviation.

## 8. Performance measures

### 8.1 Estimation

- `||beta_hat-beta_target||_1`;
- `||beta_hat-beta_target||_2`;
- maximum active-coordinate error;
- active and null coordinate bias;
- profile objective gap where an oracle optimum is available.

### 8.2 Prediction

Generate new clusters from the same DGP. Report:

- held-out pinball/check loss;
- skill relative to intercept-only and pooled QR baselines;
- quantile calibration `P(Y <= q_hat_tau)`;
- cluster-level distribution of losses.

Prediction of a new cluster excludes unavailable true random effects. A separate conditional prediction metric may be reported when repeated observations from the same test cluster are available.

### 8.3 Selection

- true-positive rate;
- false-discovery proportion and false-positive rate;
- exact support recovery;
- selected model size;
- screening sure-inclusion probability;
- stability across subject splits.

### 8.4 Inference

For active and designated null coordinates:

- bias of `beta_tilde`;
- empirical standard deviation;
- mean and median estimated standard error;
- mean-SE/empirical-SD ratio;
- 90% and 95% coverage;
- interval length;
- null type-I error;
- signal power;
- distribution and QQ plot of studentised statistics;
- null p-value uniformity diagnostics.

### 8.5 Theory mechanism

- smoothing target bias;
- effective-Hessian max-norm error;
- Schur/identity discrepancy;
- precision-row error against population/oracle direction;
- Dantzig feasibility residual;
- scaled Bahadur remainder;
- cluster-score skewness/kurtosis;
- estimated rate slope with uncertainty.

### 8.6 Computation

- total and component runtime;
- outer iterations;
- nuisance-solver iterations;
- maximum nuisance gradient;
- KKT residual;
- Hessian condition number;
- Dantzig solver status;
- peak memory;
- warning and failure class.

## 9. Pilot gate

The pilot configurations P01--P04 are diagnostic. They are considered adequate to proceed when all of the following hold:

1. geometry tests pass independently of the pilot;
2. proposed-estimator convergence rate is at least 98% in each pilot setting;
3. no systematic sign/scaling error is visible in bias or standard-error ratios;
4. median profile-identity error is below `1e-8` and maximum below `1e-6`;
5. at least 95% of requested Dantzig rows are feasible on the frozen grid;
6. mean SE/empirical SD lies roughly in `[0.80,1.20]` for the baseline, with deviations investigated rather than calibrated away;
7. 95% coverage is statistically compatible with a broad diagnostic interval such as `[0.88,0.98]` at `B=200`, accounting for MCSE;
8. increasing `n` in a targeted pilot reduces bias and Bahadur remainder;
9. no failed replication is silently dropped.

These are debugging thresholds, not final performance claims. Failing a threshold triggers a written diagnosis and code/theory review.

## 10. Final-run freeze

Before final simulations:

- tag the implementation commit;
- freeze `config/simulation_main.csv`;
- save package versions and session information;
- record hardware and worker count;
- run all geometry/unit tests;
- compute and save target approximations;
- verify that every benchmark adapter either runs or is explicitly excluded with a manuscript justification;
- generate a dry-run output and validate it against `docs/RESULTS_CONTRACT.md`.

After the freeze, changes require a new versioned run directory. Never overwrite an earlier run.

## 11. Claims permitted from simulations

The paper may claim an advantage only when supported by paired differences with Monte Carlo uncertainty and when the comparator target is commensurate. Expected patterns are:

- profile adjustment should improve inference when cluster heterogeneity is material;
- pooled modern debiasing may work when cluster effects are negligible but deteriorate as dependence strengthens;
- marginal QGEE methods may be competitive for population-average targets but answer a different question from conditional profile effects;
- oracle gaps quantify the cost of selection, nuisance estimation and precision estimation;
- under dense inverse information or severe misspecification, the method may lose coverage or produce long intervals.

Uniform dominance is neither expected nor required.
