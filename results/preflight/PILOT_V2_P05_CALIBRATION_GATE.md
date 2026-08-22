# P05 calibration gate (n=200) — verdict and round-3 evidence

Run: `results/raw/simulation/pilot_v2r2_P05b` (P05, n=200, p=500, s=5, signal=0.75,
tau=0.5, B=200/coordinate; hard calibration gate of round-2 action 14).
Pipeline: round-2 frozen implementation (unsmoothed meat primary, 2/4-fold
defect-loss mu CV, A_k/D_k diagnostics, variance ladder), commit d8ad8f5.

## 1. Gate verdict: FAIL (decisive)

| coordinate | bias | SE/SD (gate 0.80-1.20) | coverage (gate 0.88-0.98) | A_k med |
|---|---|---|---|---|
| x00001 (active) | -0.393 | 0.457 | 0.020 | 9.29 |
| x00002 (active) | +0.641 | 0.303 | 0.000 | 10.74 |
| x00006 (null)   | +0.158 | 0.730 | 0.065 | 1.38 |
| x00007 (null)   | +0.077 | 0.826 | 0.560 | 0.74 |

All coordinates feasible (1.000); mu selector concentrated at c_mu=1 (691/800),
c_mu=0.5 (108/800), c_mu=2 (1/800). The 4-fold defect-loss CV is mechanically
working (U-shaped loss, no infeasibility), but calibration is far outside the
band. Per round-2 action 14: **stop for theory review; no empirical correction
has been or will be applied** (no SE multiplier, no threshold widening, no
coverage-based selection).

## 2. Root cause: the first-stage penalised fit is identically zero

In the median replication of BOTH P01 (n=80) and P05 (n=200), the full-model
first-stage profile lasso returns beta_hat = 0 on every coordinate:

- l1_error med = 3.750 at n=80 AND n=200 (= s x signal = 5 x 0.75 = ||beta0||_1
  exactly; beta_hat is the zero vector); q25 drops (2.75 -> 1.79) but the median
  is stuck at the empty fit. max_active_error med 0.75 -> 0.487 (the active set
  is only reached in the upper tail).
- beta_hat med = 0.0000 for all four coordinates (both cells); the one-step
  beta_tilde then carries all the burden from a zero start
  (x00001: 0.35 vs target 0.75; x00002: -0.06 vs -0.75; nulls pushed to
  0.08-0.15).

Verified mechanism (standalone check, P05 config, one dataset):
- max|profile score(0)| = 0.1071 < frozen lambda = 1 x sqrt(log p / n) = 0.1763.
  The zero vector satisfies the KKT exactly (fit at multiplier 1 converges in 1
  iteration with KKT = 0).
- Frozen lambda grid (build_tuning_for_analysis_v2: multiplier fixed at 1, no
  data-driven selection): mult 0.25 -> l1_error 2.14 (nz=19, NOT converged in
  200 iters); mult 0.5 -> 3.56 (nz=3, converged); mult 1 -> 3.75 (nz=0);
  mult 2 -> 3.75 (nz=0). **No grid point reaches the assumed first-stage rate**
  (theoretical s*sqrt(log p/n) ~ 0.9).
- The score's per-coordinate sd at beta=0 is 0.022 (max 0.107 over p=500); the
  frozen lambda rate sqrt(log p/n) = 0.176 is ~8x that sd, i.e. the effective
  constant of the lambda rate is ~2-4x larger than the score scale of this DGP
  supports. The score normalisation is the theory-mandated m_i^{-1} profile
  score (AGENTS.md 5.2); nothing was changed.

## 3. The measured causal chain (n=80 -> n=200)

1. beta_hat == 0 (median rep) at both n  ->  Delta = beta_hat - beta* = -beta*,
   O_p(1), non-shrinking (l1_error 3.750 -> 3.750).
2. Dantzig residual delta = ||H omega_hat - e_k||_inf at the frozen mu anchor:
   0.284 -> 0.432 (constant order; sqrt(n)*mu_anchor =
   sqrt(log p / h) + o(1) ~ 5-6 is constant by construction).
3. Debiasing bias (e_k - H omega_hat)'Delta: |bias|/sigma0 = A_k/sqrt(n) =
   0.3015 (n=80) -> 0.2951 (n=200): **constant absolute bias, O_p(1)**.
   Per-coordinate: x00001 0.674->0.657, x00002 0.724->0.760, x00006
   0.075->0.098, x00007 0.047->0.053.
4. Hence A_k = sqrt(n)*|bias|/sigma0 grows like sqrt(n): med 2.70 -> 4.17,
   q90 6.94 -> 10.97; D_k med 11.4 -> 24.5. The Wald pivot is biased by
   ~3-11 standard errors; coverage collapses.
5. The mu-CV selector cannot rescue this: it minimises the Dantzig defect loss
   on held-out Hessians, which is blind to the zero first stage.

## 4. Oracle contrast: the variance machinery is nearly calibrated at n=200

PROFILE-DQR-TRUE-SUPPORT (lambda=0 refit on the true support):
- bias -0.008 .. +0.003; SE/SD 0.801-0.950; coverage 0.870-0.920 (two
  coordinates marginally below the 0.88 lower band).
- beta_hat med = 0.7514 (target 0.75); beta_tilde = beta_hat (one-step no-op,
  as expected at the oracle support).
- Scaling: oracle SE/SD improved 0.66-0.75 (n=80) -> 0.80-0.95 (n=200).

=> The cluster-sandwich/one-step layer is essentially closed at n=200 once the
first stage is non-degenerate. The earlier n=80 mechanism finding (large
remainder R, coarse omega) is superseded as the PRIMARY driver by the empty
first stage; a secondary ~10-20% sandwich shortfall persists even at the oracle
support (SE/SD up to 0.95, coverage up to 0.92).

## 5. Status of the other gates

- Convergence: 1.00 (median) but 15-20% of reps flag non-convergence at both n.
- Dantzig feasibility: 1.000 (round-2 selector works).
- Identity error: 0 (stable).
- P06 (n=400) still running; expected to reproduce the same empty-first-stage
  chain (lambda rate is n-independent in constant).

## 6. Questions for theory review (round-3)

Q1. The frozen first-stage rule lambda = 1 x sqrt(log p / n) (multiplier fixed,
no IC/CV) exceeds the profile-score noise floor (max|score(0)| = 0.107 vs
lambda = 0.176 at n=200; ratio ~2-4x including grid extremes), so the first
stage is identically zero in the median replication. Is the first-stage lambda
rate/constant (or its selection rule) to be revised, and if so, to what
specification (e.g., a constant relative to the empirical score scale, or a
data-driven rule)? We will implement exactly what is specified; we will not
choose lambda by coverage or truth.

Q2. The frozen mu anchor sqrt(log p/(n h)) + h^2 keeps sqrt(n)*mu constant
(~5-6) by construction, so the Dantzig slack bias is O_p(n^{-1/2}) with a large
constant. Given the round-2 amendment froze this anchor, is a first-stage fix
alone sufficient, or does the mu scale need a theory-sanctioned rate change as
well?

Q3. The oracle (TRUE-SUPPORT) coverage is 0.87-0.92 at n=200 (2/4 coordinates
below the 0.88 band) with SE/SD 0.80-0.95. Is this within the theory's
finite-sample allowance for the primary inference, or should the manuscript
state a specific higher-order correction for the cluster-sandwich at bounded
m_i (to be implemented exactly as specified)?

No code, thresholds, or tuning rules were changed for this report.
