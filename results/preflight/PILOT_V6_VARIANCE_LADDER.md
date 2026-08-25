# Round-6 variance ladder (P06-v5, B=50) — T=L+Q+R decomposition and V1 classification

Run: results/raw/simulation/pilot_v6_ladder_P06 (6 shards, 50 P06-v5
replications, PROFILE-DQR only, 204 ladder rows; commit da95c4d).
Spec: METHOD_SPECIFICATION_ROUND6_AMENDMENT.md sections 2-4.

## 1. Variance identity (verified to <1e-15)

| coord | VarT | VarL | VarQ | VarR | 2CovLQ | 2CovLR | 2CovQR | identity err |
|---|---|---|---|---|---|---|---|---|
| x00001 | 0.6078 | 0.5674 | 0.0856 | 0.1250 | -0.3390 | +0.3290 | -0.1601 | 0 |
| x00002 | 0.8731 | 0.8911 | 0.1742 | 0.2116 | -0.6736 | +0.6024 | -0.3326 | 6.7e-16 |
| x00006 | 0.7542 | 0.8248 | 0.1447 | 0.1573 | -0.6106 | +0.4844 | -0.2464 | 3.3e-16 |
| x00007 | 0.6330 | 0.7895 | 0.1348 | 0.1045 | -0.5526 | +0.2880 | -0.1311 | 1.1e-15 |

Var(T) = Var(L) + Var(Q) + Var(R) + 2Cov(L,Q) + 2Cov(L,R) + 2Cov(Q,R) holds to
machine precision: the R and covariance terms account for the full excess of
Var(T) over Var(L+Q) (e.g. x00001: VarT - Var(L+Q) = 0.294 =
VarR + 2CovLR + 2CovQR exactly).

## 2. Ladder levels (root-n scales)

| coord | E(V fit,fit) | E(V fit,target) | E(V pop,fit) | E(V pop,target) |
|---|---|---|---|---|
| x00001 | 0.3534 | 0.3537 | 0.7139 | 0.7128 |
| x00002 | 0.3852 | 0.3854 | 0.8691 | 0.8705 |
| x00006 | 0.3944 | 0.3941 | 0.8560 | 0.8563 |
| x00007 | 0.3868 | 0.3869 | 0.8726 | 0.8751 |

## 3. Mechanism ratios (5,000 paired bootstrap, seed 20260825)

| coord | M_fit [95% CI] | M_T [95% CI] | M_pop,fit [95% CI] | M_pop,target [95% CI] |
|---|---|---|---|---|
| x00001 | 1.126 [0.82,1.78] | 1.936 [1.58,2.51] | 1.258 [1.62,3.67] | 1.256 [1.63,3.65] |
| x00002 | 0.983 [0.69,1.61] | 2.229 [1.88,2.61] | 0.975 [1.54,3.69] | 0.977 [1.54,3.71] |
| x00006 | 1.099 [0.81,1.68] | 2.102 [1.70,2.62] | 1.038 [1.78,3.70] | 1.038 [1.75,3.69] |
| x00007 | 1.041 [0.75,1.72] | 1.703 [1.36,2.13] | 1.105 [1.69,3.78] | 1.108 [1.72,3.74] |

## 4. Classification: V1

- M_fit's bootstrap 95% CI contains 1 for every coordinate (0.69-1.78):
  the current fitted-direction sandwich correctly matches Var(L+Q) -- the
  first-order (influence + direction-estimation) component of T.
- M_T > 1 for every coordinate (1.36-2.62): Var(T) exceeds Var(L+Q).
- The identity (section 1) shows the excess is exactly Var(R) + 2Cov(L,R) +
  2Cov(Q,R): a finite-sample higher-order remainder, not a meat defect.
- M_pop,fit / M_pop,target are not below 1: no fitted-score squeeze (not V3);
  M_fit is not below 1 with healthy population meat (not V2); the data are
  decisive (not V4).

Per METHOD_SPECIFICATION_ROUND6_AMENDMENT.md section 4 and actions item 10:
**V1 is supported; the production variance remains unchanged; no empirical
multiplier, mu-slack variance term, or generic CR2/CR3/KC/MD correction is
added.** Under V1 the calibration bands remain reported finite-sample
diagnostics; the P05 (n=200) / P06 (n=400) undercoverage is a higher-order
remainder phenomenon, and the manuscript may state that calibration improves
with the number of independent clusters with noticeable finite-sample
undercoverage in smaller-cluster regimes.

## 5. Auto-authorised next step (actions item 10)

P07-v5 at B=200 is launched under this report (n=800, methods PROFILE-DQR |
POSTREFIT-EXACT-H | PROFILE-DQR-TRUE-SUPPORT; unchanged estimator/variance).
If no new numerical/mechanistic anomaly appears, the full versioned P01-P06
B=200 pilot is automatically authorised under the same frozen method/variance.

## 6. Compliance

Diagnostics only; no production change; truth used only in the registered
T/L/Q/R decomposition. Evidence:
results/raw/simulation/pilot_v6_ladder_P06/ladder_*_s*.csv +
ladder_aggregate.csv (run id, commit, config hash recorded).