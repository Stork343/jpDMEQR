# Pilot v2 failure decomposition (P01, corrected pipeline)

Date: 2026-08-21
Pipeline: commit 21bb8c0 (theory decisions implemented: unsmoothed primary meat,
smoothed diagnostic meat, h = c_h n^{-3/10}, Dantzig grid {0.10,...,4} with
4-fold inverse-Hessian CV + one-SE/largest-in-band selection, E1/E2/cosine/D_k/
Bahadur diagnostics, dependency-hash asset reuse).
Pilot run: results/raw/simulation/pilot_v2_P01b (P01, B=200, dev settings per design).
Gate: results/preflight/pilot_gate.json is NOT written yet (full 4-row pilot pending);
this document decomposes the P01 evidence.

## 1. What improved (vs the v1 pilot)
| Metric (P01 / PROFILE-DQR-TRUE-SUPPORT, oracle p=5) | v1 | v2 | frozen band |
|---|---|---|---|
| mean(SE)/empirical SD | 0.53-0.59 | 0.69-0.75 | [0.80, 1.20] |
| coverage (95% nominal) | 0.69-0.75 | 0.79-0.85 | [0.88, 0.98] |
| bias | ~0.004-0.008 | ~0.004-0.008 | - |

The unsmoothed quantile-score meat and the corrected bandwidth recover part of
the deficit but leave a systematic ~27-31% SE shortfall at EXACT precision rows.

## 2. Remaining deficit A: variance shortfall at the oracle support
Even with the population-correct (oracle) direction, SE/SD ~= 0.72-0.75.
The meat uses the unsmoothed score at the FITTED profile residuals:
g_i^{(0)} = -m_i^{-1} X_i' psi_tau(r_hat_i).
Hypotheses (for theory adjudication; no change made):
A1. Finite-h curvature/second-order term of the one-step expansion. With the
    new h = n^{-3/10} = 0.268 (n=80) the expected O_p(1/(nh))-style second-order
    contribution is not the right object (theory already ruled it out as a
    variance term), but the EMPIRICAL Bahadur remainder is recordable.
A2. The fitted-residual meat is "squeezed": the fitted profile residuals have
    smaller variability than the true residuals, so Var(psi_tau(r_hat)) <
    Var(psi_tau(r*)) at finite n; the sandwich understates the estimator's
    sampling variability. This is the analogue of the classic finite-sample
    sandwich under-estimation and is NOT removed by unsmoothing.
Recorded Bahadur remainder (scaled, P01): median up to ~11 in magnitude for
active coordinates -> the one-step remainder is far from negligible at n=80.

## 3. Remaining deficit B: precision rows are not first-order (D_k gate)
P01 PROFILE-DQR (full p=200):
- D_k = sqrt(n) delta_k ||beta_hat - beta_star||_1 / sigma0_pop:
  median 10.1-11.4, Q90 20.3-22.8 per coordinate (gate: median<0.5, Q90<1).
- Relative l2 row error E2_pop: median 0.55, Q90 0.72.
- cosine(POP-H direction): 0.85-0.99.
- Selected mu: c=0.5 (mu=0.284) in 614/800 rows, c=1.0 (mu=0.569) in 186/800.
- One-step correction (mean max |beta_tilde - beta_hat|): ~0.49; observed bias
  -0.30/+0.63 on the two active coordinates.

Interpretation: the selected rows satisfy feasibility (constraint residual <=
mu+1e-6) but are far from the population inverse row, so the debiasing error
(mu ||Delta||_1 ~ 0.284 * 1.4 ~ 0.4) matches the observed bias scale exactly.
The inverse-Hessian CV with the one-SE/largest-in-band rule converges to the
largest feasible multiplier (more regularised), which is precisely the value
that keeps the rows coarse. D_k quantifies the failure (median 10+, i.e. the
precision approximation is NOT in the first-order regime at n=80, p=200).

## 4. Questions for the theory side (no thresholds relaxed)
Q1 (variance): Should the primary SE use a finite-sample correction that the
   current unsmoothed meat lacks? Candidates to adjudicate (principled, not
   fitted): (a) a residual-based meat evaluated at the OBSERVED r_i without the
   fitted-squeeze (already unsmoothed); (b) an explicit small-sample term in
   the meat; (c) accept the ~0.72 ratio and document n=80 as a boundary. The
   median D_k gate already forces a fail; the question is whether ANY variance
   formula can be first-order under the present D_k regime.
Q2 (precision): Is the one-SE/largest-in-band selection the right rule at
   n=80-120? The observed selection always lands on c >= 0.5. Options to
   adjudicate: (a) smallest-loss selection (c=0.10) instead of largest-in-band;
   (b) a smaller multiplier set (0.02,0.05,0.10,...); (c) a direct row-accuracy
   target (e.g., require E2 < 0.2 by grid search) WITHOUT truth-based selection.
Q3 (oracle variance deficit): is the ~27% shortfall at exact omega attributable
   to the fitted-residual squeeze (A2) or to a term in the influence expansion
   that the current sandwich omits?

## 5. Evidence
- results/raw/simulation/pilot_v2_P01b/coordinate_metrics.csv (per-rep D_k,
  bahadur_scaled, E1/E2/cosine, mu, se/estimated_se_smoothed)
- results/raw/simulation/pilot_v2_P01b/replication_metrics.csv
- results/tables/pilot_v2_P01b/coverage_summary.csv and theory_summary.csv
- results/intermediate/population_directions/P01.rds (h_analysis=0.2686,
  sigma0_pop stored)
