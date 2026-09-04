# P05 round-3 lambda-mechanism check (B=50) — verdict: STOP, escalate

Run: `results/raw/simulation/pilot_v3_P05_B50` (P05, n=200, p=500, s=5,
signal=0.75, tau=0.5; B=50, three shards merged; commit 50fedbe).
Spec: `docs/METHOD_SPECIFICATION_ROUND3_AMENDMENT.md` (cluster self-normalised
profile-score loadings, two-pass, weighted normalised KKT, solver budgets
2000/50/1e-7). Decision gate: `results/preflight/PILOT_GATE_THEORY_ACTIONS.txt`
items 13-14.

## 1. Verdict

The round-3 lambda rule is implemented exactly as specified and works
mechanically, but the one-step diagnostics are **not coherent**:
- PROFILE-DQR active-coordinate coverage 0.000 / 0.020 (x00001 / x00002),
  bias -0.335 / +0.508, A_k 6.7 / 9.4, scaled Bahadur remainder -13.1 / +17.8.
- Null coordinates are near calibration (coverage 0.36 / 0.88, A_k 0.74 / 0.48),
  which isolates the failure to the active-coordinate first-stage shrinkage.
- Per action 14, the full P01-P06 rerun is **NOT started**; this is a STOP-and-
  escalate record. No code, threshold, or tuning rule was changed.

## 2. What the round-3 fix achieved (mechanical pass)

- First stage non-degenerate: l1_error med 3.266 (q25 3.186, q75 3.344) vs the
  old empty fit 3.750; first_stage_nonzero_count med 3; tpr med 0.60.
- KKT acceptance 100%: final_kkt_normalized med 1.3e-4, p90 7.1e-4 (<= 1e-3);
  status ok = 1.000 for all 50 replications; iterations med 23 (budget 2000).
- Zero fit is infeasible: zero_kkt_ratio med 1.91 (was <1 under the old rule).
- lambda_rule=cluster-score-loading-v3 audit fields populated in every row;
  loading pass0/pass1 medians stable (0.2143 / 0.2137), lambda_base=0.3205,
  median coordinate penalty 0.0692.

## 3. The remaining gap (sharpened mechanism)

The score-calibrated penalty sits at the **identification boundary** for the
active coordinates:

- Median max|profile score(0)| = 0.128 vs lambda_0 = 0.3205 (zero infeasible by
  ~2x), but the ACTIVE-coordinate score at beta=0 is only ~1.0-1.5x its
  coordinate penalty (median coordinate penalty 0.069 vs signal 0.75):
  the profile effective curvature in this DGP (AR1 rho_x=0.5, h=n^{-3/10}=0.204,
  bounded m_i, q=1 intercept nuisance) is weak enough that the calibrated
  penalty dominates the signal at n=200.
- Result: the median replication zeroes x00002 entirely and shrinks x00001 to
  ~0.22 (target 0.75); Delta = beta_hat - beta* keeps an O_p(1)-scale active
  component, and the one-step inherits it (A_k 6.7-9.4, i.e. bias ~7-9 standard
  errors). The observed l1_error 3.27 is ~3.6x the retained rate
  s*sqrt(log p/n) ~= 0.9.
- POP-H (population direction, same first stage) halves the active bias
  (coverage 0.54 / 0.34, A_k 4.2 / 5.2) but remains out of band; TRUE-SUPPORT
  (oracle refit) is near calibration (coverage 0.86-0.92), confirming the
  variance layer is not the failing part.

## 4. Questions for theory review (round-4)

Q1. The retained first-stage rate assumes profile restricted strong convexity
with the score-calibrated lambda. Observed: the active-coordinate score at
beta=0 is comparable to the calibrated penalty at n=200, so the first stage
selects only ~60% of the active set with heavy shrinkage. Is the curvature
regime of this DGP (AR1 rho_x=0.5, h=c_h n^{-3/10}, bounded m_i, q=1) at
n=200 consistent with the rate, and at what n does identification engage?

Q2. The primary multiplier is frozen at a=1. The active-score-to-penalty ratio
~1 suggests a theory-sanctioned revision of a (or a recorded data-driven
selection along the registered a-grid) may be required. We will implement
exactly what is specified; we will not select a by coverage or truth.

Q3. Is the analysis bandwidth h=n^{-3/10} compatible with first-stage
identification at n=200? h-smoothing attenuates the score/curvature scale.

Q4. Optional cheap scaling probe: run P06 (n=400) at B=50 (~7-10h on 3
shards, machine currently idle) to measure whether l1_error shrinks
(3.27 -> ?) and A_k drops (6.7 -> ?) before committing ~2 days to the full
pilot. Authorise?

## 5. Compliance

No empirical multiplier, threshold widening, truth/coverage-based selection,
or score/Hessian/target change. The round-3 specification was implemented
verbatim; the remaining gap is an identification/curvature question for the
theory side.

Evidence: results/raw/simulation/pilot_v3_P05_B50/{replication_metrics,
coordinate_metrics,theory_diagnostics,screening_records}.csv + merge_manifest.json.
