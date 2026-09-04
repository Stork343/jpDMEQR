# Round-4 paired mechanism probe (P05-v4 / P06-v4, B=50) — report

Runs: results/raw/simulation/pilot_v4_P05_B50 (n=200) and pilot_v4_P06_B50
(n=400), merged from 3 shards each, 0 runner errors, nested clusters (P06 uses
P05 training/test seeds; seed_experiment=P05 recorded in run_request.json).
Commit a446a76. Spec: METHOD_SPECIFICATION_ROUND4_AMENDMENT.md.
Per PILOT_GATE_THEORY_ACTIONS.txt item 15: diagnostic only, no formal gate.

## 1. First stage: the dual-bandwidth fix engages as designed

| metric (PROFILE-DQR, B=50) | P05-v4 n=200 | P06-v4 n=400 |
|---|---|---|
| l1_error med (q25, q75) | 3.243 (3.139, 3.312) | 2.273 (2.183, 2.388) |
| paired l1 diff (nested) | — | -0.957 |
| nonzero count med / tpr med | 3 / 0.60 | 5 / 1.00 |
| exact support | 0.040 | 1.000 |
| KKT normalized med / p90 | 1.2e-4 / 4.6e-4 | 1.7e-4 / 6.6e-4 |
| status ok | 1.000 | 1.000 |
| h_est / h_est_over_rsc | 0.4199 / 1.065 | 0.3531 / 1.267 |
| coord penalty med | 0.0672 | 0.0485 |

RSC diagnostics healthy and improving: lambda_min(SS) fitted 0.060 -> 0.079,
target 0.105 -> 0.110, population 0.110 -> 0.112, cone proxy 0.194 -> 0.220,
condition 10.1 -> 8.0. At n=400 the first stage converges with the exact
support in every replication. The n=200 cell remains at the identification
boundary (tpr 0.60), exactly as the round-4 order-level analysis predicted.

## 2. Debiased layer: improving, but far from the gate at both n

PROFILE-DQR (active coords):
- x00001: A_k 6.53 -> 4.19; bias med -0.347 -> -0.143; coverage 0.040 -> 0.280
  (mcse 0.028 -> 0.063); Bahadur -13.0 -> -10.9.
- x00002: A_k 9.40 -> 8.35; bias med +0.561 -> +0.329; coverage 0.040 -> 0.060
  (mcse 0.028 -> 0.034); Bahadur +17.7 -> +17.2.
- Nulls: coverage 0.44/0.86 -> 0.70/0.96.

POP-H: x00001 cov 0.540 -> 0.860, x00002 cov 0.340 -> 0.700; A_k 4.2/5.2 ->
2.0/3.2. TRUE-SUPPORT: 0.86-0.92 -> 0.88-0.96.

The trajectory is correct and the machinery is coherent, but the practical
P05/P06 calibration is not achieved at either n.

## 3. Sharpened quantitative observation: the dominant remaining bias is the
frozen mu-slack term, and it decays too slowly

At n=400 the first stage has the EXACT support, yet x00002 still carries
bias +0.329 (A_k 8.35). Decomposing the one-step bias:

- A_k/sqrt(n) = |(e_k - H omega_hat)' Delta|/sigma0: x00002 0.67 (n=200) ->
  0.42 (n=400); x00001 0.46 -> 0.21. The absolute bias shrinks only like
  n^{-1/2}, so the bias-to-SE ratio A_k stays order 1 (slowly decreasing).
- The Dantzig slack contribution is bounded by sqrt(n)*mu*||Delta||_1/sigma0
  with sqrt(n)*mu = sqrt(log p / h_inf) + o(1) ~ 6-7 (constant by the frozen
  anchor). With ||Delta||_1 = O_p(s sqrt(log p/n)) this contributes
  A_k_slack ~ s log p / (sigma0 sqrt(n h_inf)): ~4.9 (n=200), ~3.8 (n=400),
  ~3.0 (n=800), ~1.5 (n=5400) — an n^{-0.35}-type decay that keeps the bias
  comparable to (or above) the parametric SE at every practical n.

This matches the observed A_k trajectory (6.5->4.2, 9.4->8.4) and the
persistent Bahadur remainders (+-11 to +-18). Round-3 Q2 flagged exactly this
term; round-4 kept the mu anchor unchanged, with the reopen condition "if,
after lambda repair, A_k, POP-H row errors, total Bahadur remainder and
P05-P06 scaling remain poor jointly". The probe shows: first-stage/RSC fixed,
POP-H row errors much improved, but A_k and Bahadur remain poor with slow
n-scaling.

## 4. Recommendation (for theory/project review, round-5)

Per action item 16, the full P01-P06 rerun is NOT launched on the basis of
this report. Questions for adjudication:

Q1. The mu-slack bias term (sqrt(n)*mu constant ~6-7 by the frozen anchor)
appears to be the dominant residual one-step bias at exact support, with an
n^{-0.35} decay that cannot bring coverage into band at n=200/400 (nor
realistically at n=800). Is the round-4 reopen condition met, and if so what
mu-rate amendment is sanctioned (exact formula)?

Q2. If Q1 is deferred: authorise P07 (n=800, B=200) directly as the scaling
locator without first running the full P01-P06, given the probe's consistent
evidence that the n=200 practical gate cannot pass? (The Q5 conditional path
was written as "after the full P05 gate still fails"; the B=50 evidence makes
the full-run outcome near-certain.)

Q3. The one-step at exact support still under-corrects (x00002 bias +0.33 at
n=400, se/sd 0.26). Is the residual within the theory's leading order given
the mu term, or does the h_inf-based direction require a theory-sanctioned
change?

## 5. Compliance

No code/threshold/tuning change since a446a76. All round-3/4 rules implemented
verbatim. All 300 probe replications (50 x 2 cells x 5 methods... 250 method
rows per cell) accounted for, 0 runner errors, KKT 100%.

Evidence: results/raw/simulation/pilot_v4_P05_B50/{replication_metrics,
coordinate_metrics,theory_diagnostics,screening_records}.csv +
results/raw/simulation/pilot_v4_P06_B50/ (both with merge_manifest.json).
