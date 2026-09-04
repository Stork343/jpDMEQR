# Pilot v5 final report — P01-P06 B=200 under the frozen round-5/6 method

Run: results/raw/simulation/pilot_v5_B200 (6 shards merged; 7600 method rows,
24000 coordinate rows, 1200 screening records; 0 runner errors; commit
1dcc159-era code with portable launcher).
Method (frozen): score-loaded L1 at h_est=(log p_P/n)^{1/4} -> post-L1 profile
refit at h_inf=n^{-3/10} -> full-p reprofile -> Dantzig (round-2 CV rule) ->
one-step -> unsmoothed fitted cluster-score sandwich. Variance: round-6 V1
classification (no production change).

## 1. Overall

All 6 cells: status ok = 1.000, converged = 1.000, no failure stages, KKT
acceptance and post-refit acceptance 1.000. No implementation anomaly.

## 2. Per-cell PROFILE-DQR results (B=200, coverage mcse +/-0.02-0.03)

| cell | n | error/tau/q | l1_refit | active sel | coverage range | SE/SD range |
|---|---|---|---|---|---|---|
| P01 | 80 | gaussian | 2.17 | 3/5 | 0.26-0.68 | 0.41-0.55 |
| P02 | 120 | t3, q=2 | 4.38 | 2/8 | 0.15-0.75 | 0.37-0.60 |
| P03 | 120 | skew tau=.25 | 4.45 | 2/8 | 0.15-0.73 | 0.45-0.57 |
| P04 | 120 | contam tau=.75, q=2 | 4.61 | 2/8 | 0.21-0.68 | 0.41-0.52 |
| P05 | 200 | gaussian (gate) | 1.59 | 3/5 | 0.37-0.75 | 0.27-0.62 |
| P06 | 400 | gaussian | 0.23 | 5/5 | 0.79-0.885 | 0.66-0.78 |

P07-v5 B=200 (n=800): l1_refit 0.158, coverage 0.865-0.910 (all coverage CIs
contain the lower band edge 0.88), SE/SD 0.75-0.83. TRUE-SUPPORT at n=800:
0.92-0.935 (in band).

## 3. Interpretation (round-6 V1 framework)

- The T=L+Q+R ladder (P06, B=50) verified sand-which == Var(L+Q) (M_fit CI
  contains 1 for all coords) and Var(T)-Var(L+Q) = Var(R)+2Cov(L,R)+2Cov(Q,R)
  to machine precision: the missing coverage is a finite-sample higher-order
  remainder, not a first-order defect.
- Coverage improves monotonically with n: 0.15-0.75 (n=120) ->
  0.37-0.75 (n=200) -> 0.79-0.885 (n=400) -> 0.865-0.91 (n=800); SE/SD
  0.27-0.68 -> 0.62-0.80 -> 0.75-0.83.
- P02-P04 (first cells to exercise t3/skewness/contamination and q=2): the
  small-n selection recovers only 2/8 actives; the post-refit cannot repair
  omitted signal (bias +0.17..0.22 on the missed active) -- the expected
  model-selection phenomenon stated in the round-5 amendment.
- P05 (n=200, the hard gate cell) fails the band at B=200 (coverage
  0.37-0.75, CI upper 0.81 < 0.88): under V1 this is reported as a
  finite-sample scaling fact, not an implementation defect.

## 4. Manuscript statements supported

- First-order asymptotic machinery: fully supported (rates, one-step
  expansion, sandwich correctness per the V1 ladder).
- Finite sample: "calibration improves with the number of independent
  clusters, with noticeable finite-sample undercoverage in the smaller-cluster
  regimes"; the n-scaling table above is the evidence.
- No nominal-coverage claim at n <= 400; no high-dimensional confirmatory Wald
  at n=129 (GSE65391, round-6 decision, now empirically documented via
  PAPP-HD: coverage 0.33-0.76 with 2/5 active-selection misses).

## 5. Freeze status

- Pilot record complete (this report + raw outputs + merge manifests).
- mechanism_evidence (formal_gate=false) list under pilot_gate.json:
  P05/P06-v4 B=50, P05/P06-v5 B=50, P07-v5 B=50, round-6 ladder, full
  P01-P06 B=200, P07-v5 B=200.
- Next: final pilot manifest + freeze-preflight here; B=500/1000 main registry
  (compute: cloud server).

## 6. Compliance

No empirical multiplier, no threshold widening, no DGP weakening, no silent
deletion of failures (all 1200 replications present), solver/estimator/
variance unchanged since the round-6 ladder classification.