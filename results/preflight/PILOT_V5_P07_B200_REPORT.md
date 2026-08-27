# P07-v5 B=200 (n=800) — V1-path scaling result

Run: results/raw/simulation/pilot_v5_P07_B200 (6 parts: s1, s3-s6, s2r catch-up
after the s2 OOM; merged 600 method rows = 200 tasks x 3 methods; 0 runner
errors; commit 637e4e1 + adapter fix commit).
Spec: METHOD_SPECIFICATION_ROUND6_AMENDMENT.md section 5 / actions item 10.

## 1. First stage

l1_lasso med 1.530 -> l1_refit med 0.158 (l2 0.072); selected_size med 5
(all actives), refit_set med 7, omitted_active 0, post_refit ok 1.000,
KKT p90 7.6e-4, status ok 1.000, no failure stages.

## 2. Practical method at n=800 (B=200)

| coord | bias | coverage (mcse) | SE/SD | band-consistent* |
|---|---|---|---|---|
| x00001 | -0.0005 | 0.880 (0.023) | 0.83 | yes |
| x00002 | -0.0026 | 0.865 (0.024) | 0.75 | yes |
| x00006 | -0.0006 | 0.870 (0.024) | 0.77 | yes |
| x00007 | -0.0001 | 0.910 (0.020) | 0.83 | yes |

*95% binomial CI on coverage contains the lower band edge 0.88 (CIs:
[0.835,0.925], [0.818,0.912], [0.823,0.917], [0.870,0.950]).

TRUE-SUPPORT: coverage 0.920-0.935, SE/SD 0.92-0.96 (in band).
POSTREFIT-EXACT-H (diagnostic; see fix below): coverage 0.95-0.96,
SE/SD 1.04-1.07.

The three-cell practical scaling (n=200/400/800) under the round-five method:
coverage 0.36-0.82 -> 0.72-0.86 -> 0.865-0.910; SE/SD 0.27-0.68 -> 0.62-0.80
-> 0.75-0.83. Monotone improvement; at n=800 the coverage sits at/just below
the band edge, statistically consistent with the lower edge, and is
attributed per the round-six ladder to the V1 higher-order remainder.

## 3. Diagnostic fix transparency

The POSTREFIT-EXACT-H 95% intervals used zcrit = qnorm(1 - c/2) = 0.063
instead of qnorm(1 - (1-c)/2) = 1.96, making those intervals ~31x too narrow
(coverage 0.05-0.065 instead of 0.95-0.96). This affected only the
simulation-only POSTREFIT-EXACT-H diagnostic; the primary PROFILE-DQR
intervals go through debias_profile_coordinates_v2 (correct zcrit). The
adapter was fixed, a regression test asserts width = 3.92*se, and the merged
B=200 coordinate metrics for POSTREFIT-EXACT-H rows were recomputed offline
from the stored (correct) se. No primary result changed.

## 4. V1-path decision (actions item 10)

No new numerical/mechanistic anomaly in the practical estimator: 200/200
replications ok, KKT and post-refit acceptance 1.000, bias <= 0.003, monotone
scaling toward the band. The POSTREFIT-EXACT-H zcrit defect was a documented
diagnostic-only bug now fixed and tested. Per action 10, the full versioned
P01-P06 B=200 pilot under the same frozen method/variance (round-5 estimator,
round-6 unchanged primary sandwich) is therefore authorised.

## 5. Compliance

Gates remain reported diagnostics under V1; no threshold, estimator, or
variance change; the diagnostic fix is documented. Evidence:
results/raw/simulation/pilot_v5_P07_B200/{replication_metrics,
coordinate_metrics,theory_diagnostics,screening_records}.csv +
merge_manifest.json (parts s1, s3, s4, s5, s6, s2r).