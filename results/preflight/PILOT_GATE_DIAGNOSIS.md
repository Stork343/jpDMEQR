# Pilot gate diagnosis (P01-P04, B=200, dev settings per design)

Date: 2026-08-19
Run: results/raw/simulation/pilot_final_merged (partitions P01-P04 merged)
Commit: dfcc379 (R code identical to the run-time 36527be; CI-only commits)
Gate manifest: results/preflight/pilot_gate.json (pass = FALSE)

## Gate outcome
- convergence_pass      TRUE
- identity_pass         TRUE  (profile identity ~1e-16 at analysis scale)
- dantzig_pass          TRUE  (feasibility >= 0.95)
- no_runner_failures    TRUE
- no_missing_methods    TRUE
- se_ratio_pass         FALSE (baseline P01: se_sd_ratio 0.50-0.87 vs [0.80, 1.20])
- coverage_pass         FALSE (baseline P01: coverage 0.005-0.86 vs [0.88, 0.98])

## Finding A (primary): sandwich variance underestimates the empirical SD by ~40-50%
Even the ORACLE-SUPPORT adapter (p = 5, no screening, no lasso shrinkage,
bias ~0.004 on all coordinates) shows se_sd_ratio ~0.53-0.59 and coverage
~0.69-0.75 at n = 80 (P01). The shortfall is therefore NOT caused by the
high-dimensional Dantzig rows, the debias correction, or screening.

Implementation of cluster_sandwich_coordinate_v2 matches METHOD_SPEC sec.5.5
exactly (centred projected cluster scores, mean square, se = sqrt(sigma2 / n),
df_correction = FALSE by default). The (n-1)/(n) correction would change the
ratio by ~1% only, so no finite-sample multiplier in the frozen formula can
explain a 1.8x SE shortfall.

Hypotheses for the theory side (no changes made):
1. Smoothing width h = n^{-1/3} makes Var(psi_h) < tau(1-tau), so the smoothed
   sandwich SE underestimates the sampling variability of the estimator
   (which is driven by the unsmoothed residual mass near the kink).
2. Finite-sample curvature of the profiled criterion (the sandwich ignores
   the estimation noise of the first-stage beta_hat and gamma_hat beyond the
   first-order expansion; the second-order term may be non-negligible at
   n = 80-120 with cluster sizes 2-8).
3. Cluster-size imbalance (m_i in 2:8): the m_i^{-1} weighting in the score
   and the sandwich may need a size-weighted finite-sample correction.

## Finding B (secondary): debias bias at the full design (n = 80, p = 200)
PROFILE-DQR at P01: bias -0.26 (x00001) / +0.63 (x00002), coverage 0.005-0.37.
mean_precision_row_error ~3.9 and the selected Dantzig mu ~0.29
(mu = c_mu * [(log p / (n h))^{1/2} + h^2], smallest feasible frozen grid
value) imply the precision row omega is far from the true inverse row, so
the one-step correction (mean one_step_correction ~0.49) is applied but
inaccurate. TRUE-SUPPORT (exact omega) removes the bias (bias ~0.004), which
confirms the bias originates in the precision rows, not in the target, the
DGP, or the score definition.

## What this means
- The pilot gate did its job: the failure was detected BEFORE the final
  registry (B = 500/1000) was launched, saving days of compute.
- The failures concern the frozen variance formula (METHOD_SPEC sec.5.5) and
  the frozen Dantzig mu grid. Both are theory-side decisions (AGENTS.md
  sec.18): do not adjust thresholds, mu, or variance formulas without
  theory review.
- The population target assets, direction assets, benchmark acceptance,
  geometry, and micro-preflight are unaffected by these findings and remain
  valid (they do not depend on the sandwich or the analysis mu grid).

## Evidence files
- results/tables/pilot_final_merged/coverage_summary.csv
- results/tables/pilot_final_merged/theory_summary.csv
- results/tables/pilot_final_merged/method_summary.csv
- results/raw/simulation/pilot_final_merged/* (raw per-rep tables)
