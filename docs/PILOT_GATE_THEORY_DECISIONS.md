# Pilot-gate theory decisions

The current authoritative second-round decision is:

`results/preflight/PILOT_V2_THEORY_DECISION_ROUND2.md`

It supersedes the Dantzig selection/gating and pilot-calibration provisions of `results/preflight/PILOT_GATE_THEORY_DECISION.md`. All non-conflicting round-one decisions remain in force.

The following rules are now authoritative:

1. The inferential target remains the **unsmoothed regularised profile parameter**. The exact profile score and Schur-complement Hessian definitions in `docs/METHOD_SPECIFICATION.md` are unchanged.
2. The primary asymptotic Wald meat remains the ordinary quantile cluster score at fitted profile residuals; smoothed-score meat remains a diagnostic. The residual oracle SE/SD shortfall is **not** repaired by a scalar multiplier. A target-score/fitted-score/population first-order variance ladder and a TRUE-SUPPORT delete-one-cluster jackknife must determine whether the remaining gap is score plug-in/leverage bias or a genuine nonlinear remainder.
3. No additional scalar correction for unequal `m_i` is authorised. P01's bounded non-informative cluster sizes are already represented by `m_i^{-1}` cluster scores.
4. Primary inference bandwidth remains `h=c_h n^{-3/10}`, `c_h in {0.75,1,1.25}`.
5. The previous **one-SE/largest-in-band Dantzig rule is revoked**. For `n<200`, use deterministic two-fold cluster CV; for `n>=200`, use four folds. The primary validation loss is held-out inverse defect `||H_val omega_train-e_k||_inf`. Select the smallest mean defect; numerical ties go to the **smaller mu**. The inverse-quadratic loss is diagnostic only.
6. The Dantzig multiplier grid for the next pilot is `c_mu in {0.02,0.05,0.10,0.25,0.50,1,2,4}`. Infeasible candidates are recorded; no further downward expansion is allowed without another theory decision.
7. `D_k=sqrt(n) delta_k ||beta_hat-beta_star||_1/sigma0_pop` is a conservative Hölder upper bound. The old `median(D_k)<0.5, Q90(D_k)<1` values are **not hard freeze gates**. Record `D_k`, but also record the actual simulation-only normalized inverse-defect inner product, POP-H row errors/cosine similarity and the normalized Bahadur remainder. None may tune `mu`.
8. P01--P04 remain stress/diagnostic cells. Add P05 (`n=200,p=500,s=5,tau=0.5,q=1`, baseline Gaussian random intercept) and P06 (`n=400,p=500,s=5`, same DGP). The `[0.80,1.20]` SE/SD and broad coverage calibration thresholds are **not relaxed**; the hard first-order calibration gate is assessed at P05, with P06 used to confirm scaling. P01--P04 must still pass correctness/identity/failure-accounting gates and are reported even when calibration is poor.
9. Profile-target assets remain reusable when their dependency hash is unchanged. Their own small target-approximation bandwidth need not equal the analysis bandwidth. POP-H assets are defined at `h_analysis=c_h n_analysis^{-3/10}`; the population size controls Monte Carlo accuracy only. Top-level `n_analysis`, `h_analysis`, `sigma0_pop`, and `Sigma0_population` are the intended audit fields.
10. HiGHS is authorised as an LP solver for the identical Dantzig program after a prospective parity suite against CLARABEL on at least 20 representative rows. The solver change is numerical, not methodological.
11. Keep Module-A `n=100` as a prespecified finite-cluster boundary cell. Manuscript language is result-contingent: calibration is assessed over the `n=100,200,400` sequence; say it approaches nominal only if the results show that trend.
12. GSE65391 remains viable, but confirmatory inference at roughly 129 subjects stays blocked until the small-sample variance mechanism is resolved. Low-dimensional confirmatory modules require a prespecified subject-level jackknife/leverage-correction audit; high-dimensional exploratory gene intervals require acceptable precision-row stability in a comparable simulation regime.

No final-scale `B=500/1000` simulation or confirmatory empirical inference is authorised until the round-two implementation/action sequence is completed and a fresh freeze manifest passes.