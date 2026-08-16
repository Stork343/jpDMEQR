# Benchmark matrix and implementation contract

The final simulation study must include methods that represent four distinct comparison classes: pooled high-dimensional quantile regression, classical mixed quantile regression, recent longitudinal high-dimensional quantile regression, and oracle versions of the proposed method. A competitor is not considered implemented until its adapter passes its own unit test on a small known dataset.

| ID | Method | Role | Target/structure | Dimension support | Inference compared? | Current repository status |
|---|---|---|---|---|---|---|
| `PROFILE-DQR` | Corrected profile estimator + row-wise Dantzig + cluster sandwich | Proposed | Regularised conditional profile target | High-dimensional after optional split screening | Yes | Reference implementation in `R/profile_v2.R`; optimization required |
| `POOLED-QR-LASSO` | Pooled l1 quantile regression | Classical negative control | Ignores clustering | High-dimensional | Estimation/prediction; naive inference only if clearly labelled | Adapter implemented |
| `SQR-DEBIASED-IID` | Debiased smoothed QR for independent observations | Modern negative control | Independent-observation target; clustering ignored | High-dimensional | Yes, labelled naive under clustering | Adapter specification present; full parity with Yan et al. required |
| `LQMM` | Geraci/Bottai linear quantile mixed model | Classical direct comparator | Conditional quantile mixed model via AL working likelihood | Low/moderate dimension | Fixed-effect intervals where stable | Adapter implemented for common/true support |
| `BIAS-ADJ-LQMM` | Bias-adjusted clustered quantile estimator | Many-small-cluster comparator | Fixed-effect target with bootstrap bias adjustment | Low/moderate dimension | Yes | `implementation_required` |
| `DOUBLE-PEN-QLMM` | Li, Liu & Luo double-penalised QLMM | Direct penalised mixed-QR comparator | Fixed and random-effect selection | Moderate/high dimension per paper | Estimation/selection; inference if reproduced faithfully | `implementation_required` |
| `QGEE-SCAD` | Zu, Lian, Green & Yu quantile penalised GEE with SCAD | Recent longitudinal comparator | Marginal/estimating-equation quantile target | Ultra-high dimensional with screening | Estimation/selection and its reported asymptotic SE if implementable | `implementation_required` |
| `QIF-SEE` | Bhattacharya, Bhuiyan & Chatla QIF + smooth-threshold equations | Recent SJS comparator | Longitudinal quantile selection target | Diverging dimension | Estimation/selection; not treated as equivalent profile inference | `implementation_required` |
| `BAYES-MIXED-LASSO` | Bayesian Gaussian mixed model with shrinkage | Supplemental prediction comparator | Conditional mean, not quantile | Moderate dimension after screening | No coverage comparison | Existing legacy adapter; must be relabelled |
| `PROFILE-DQR-TRUE-SUPPORT` | Proposed method on true support | Oracle | Same profile target | Oracle low dimension | Yes | Implemented through `penalty_factor`/column subset |
| `PROFILE-DQR-TRUE-NUISANCE` | Proposed score/curvature using generated random effects | Oracle diagnostic | Structural nuisance known | Oracle | Yes | `implementation_required` |
| `PROFILE-DQR-POP-H` | Proposed estimator with population effective direction | Oracle diagnostic | Same target | Oracle | Yes | Interface specified; target-specific implementation required |
| `ORACLE-LQMM` | LQMM on true support | Oracle classical | LQMM working target | Oracle low dimension | Yes where stable | Adapter implemented |

## Required citations/identifiers

- Geraci and Bottai (2014), *Statistics and Computing*, DOI `10.1007/s11222-013-9381-9`.
- Geraci (2014), `lqmm` package paper, DOI `10.18637/jss.v057.i13`.
- Li, Liu and Luo (2020), *Journal of Systems Science and Complexity*, DOI `10.1007/s11424-020-9065-4`.
- Battagliola et al. (2022), *Econometrics and Statistics*, DOI `10.1016/j.ecosta.2021.07.003`.
- Zu, Lian, Green and Yu (2023), *JASA*, DOI `10.1080/01621459.2022.2128806`.
- Yan, Wang and Zhang (2023), *JMLR* 24(245), paper identifier `22-1217`.
- Bhattacharya, Bhuiyan and Chatla (2026), *Scandinavian Journal of Statistics*, DOI `10.1111/sjos.70077`.

## Adapter API

Every benchmark adapter must expose the same minimum interface:

```r
fit_benchmark_<id>(
  train,
  tau,
  target_coords,
  tuning,
  seed,
  control = list()
)
```

and return a named list containing:

```text
method_id
status                 ok / warning / failed / not_implemented
beta_hat               named numeric vector on the supplied design
beta_tilde             optional named numeric vector
se                      optional named numeric vector
ci_lower, ci_upper      optional named vectors
selected                integer or character feature identifiers
prediction_function     function or serialisable prediction object
objective
runtime_sec
converged
kkt_residual
warning_messages
failure_stage
implementation_version
```

Adapters must never return a successful status after catching an error and substituting zeros or another method.

## Fair comparison by task

### Estimation and selection

Compare all methods that estimate a commensurate fixed-effect vector. Low-dimensional methods receive either the true support (oracle table) or the same frozen screened set (practical table). Distinguish these tables.

### Prediction

All fitted methods can be compared by held-out pinball loss if they produce quantile predictions for new subjects. Mean-model Bayesian competitors are reported separately and are not interpreted as quantile estimators.

### Inference

Primary coverage comparisons include methods with an explicit coordinate-wise inferential procedure. If a competitor's target is marginal while the proposed target is conditional/profiled, show both methods but state that coverage is evaluated for each method's own target or restrict the comparison to designs where targets coincide.

### Variable selection

SCAD/QIF procedures can be compared on TPR/FDP/model size. Selection accuracy is not a substitute for interval validity.

## Implementation acceptance tests

Before a benchmark enters the main registry:

1. It fits a low-dimensional Gaussian random-intercept example without error.
2. It returns coefficient names aligned with the supplied design.
3. Its prediction function accepts new clusters without using generated latent random effects unless the task is explicitly conditional prediction.
4. Its tuning uses training clusters only.
5. Its output includes solver status and runtime.
6. A deterministic seed reproduces identical output.
7. At least one limiting-case check agrees with a simpler method, such as LQMM approaching ordinary quantile regression when random-effect variance is negligible.

## Missing implementation policy

The config files may list an adapter whose status is `implementation_required`. The runner must stop before a final run if any required method is missing. During a pilot, `--allow-missing-benchmarks=true` may be used, but the output must record the omission and no final comparative claim may be produced.
