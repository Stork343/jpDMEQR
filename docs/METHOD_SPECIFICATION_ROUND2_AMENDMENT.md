# Method specification — round-two pilot amendment

**Status:** authoritative amendment to `docs/METHOD_SPECIFICATION.md` pending incorporation into the base file  
**Date:** 2026-08-21  
**Scope:** row-wise precision tuning, precision diagnostics, small-sample variance diagnostics, and pilot-calibration regime only  
**Source decision:** `results/preflight/PILOT_V2_THEORY_DECISION_ROUND2.md`

This amendment supersedes the conflicting sentences in Sections 6 and 9.4 of `docs/METHOD_SPECIFICATION.md` that refer to (i) selecting the smallest feasible Dantzig value without validation, or (ii) selecting the largest value within one standard error of an inverse-quadratic CV minimum. The Dantzig linear program itself, exact profile score/Hessian, one-step estimator, primary unsmoothed meat, and `h=c_h n^{-3/10}` remain unchanged.

## 1. Dantzig candidate grid

For target coordinate `k`, solve

\[
\widehat\omega_k(\mu)
\in\arg\min_{\omega}\|\omega\|_1
\quad\text{subject to}\quad
\|\widehat H_{\mathrm{eff}}\omega-e_k\|_\infty\le\mu.
\]

The candidate tolerance is

\[
\mu(c)=c\left\{\sqrt{\log(p)/(nh)}+h^2\right\},
\qquad h=c_h n^{-3/10},
\]

with the frozen second-round grid

\[
c\in\{0.02,0.05,0.10,0.25,0.50,1,2,4\}.
\]

Infeasible candidates are recorded. The grid may not be expanded below `0.02` after viewing coverage, bias, or truth-based row errors.

## 2. Non-truth selection of `mu`

Use deterministic cluster-level cross-validation:

- `n < 200`: two folds;
- `n >= 200`: four folds.

For each fold and candidate, solve the row on the training-cluster Hessian and evaluate the held-out inverse defect

\[
L_{f,k}^{\mathrm{defect}}(c)
=\|H_f\widehat\omega_{-f,k}(c)-e_k\|_\infty.
\]

The primary selected candidate minimises the mean held-out inverse defect. If candidate losses are tied within a committed numerical tolerance, choose the **smaller `mu`**. Do not use a one-standard-error enlargement. Re-solve the selected row on the full-sample Hessian.

The inverse-quadratic loss

\[
\frac12\omega^\mathsf T H_f\omega-e_k^\mathsf T\omega
\]

is retained as a secondary diagnostic only.

No simulation truth, POP-H direction, empirical bias, interval coverage, or p-value may enter the selection rule.

## 3. Precision diagnostics

Continue to report the Dantzig feasibility residual and row norms. Where a population direction is available, additionally report

\[
E_{2,k}=\frac{\|\widehat\omega_k-\omega_k^{\mathrm{pop}}\|_2}
{\|\omega_k^{\mathrm{pop}}\|_2},
\qquad
E_{1,k}=\frac{\|\widehat\omega_k-\omega_k^{\mathrm{pop}}\|_1}
{1+\|\omega_k^{\mathrm{pop}}\|_1},
\]

cosine similarity, and

\[
D_k=\frac{\sqrt n\,\delta_k\|\widehat\beta-\beta^\star\|_1}
{\sigma_{0,k}^{\mathrm{pop}}}.
\]

`D_k` is a conservative Hölder upper bound and is **not a hard freeze threshold**. Also report the actual simulation-only inverse-defect contribution

\[
A_k=\frac{\sqrt n\,
|(e_k-\widehat H\widehat\omega_k)^\mathsf T
(\widehat\beta-\beta^\star)|}
{\sigma_{0,k}^{\mathrm{pop}}}.
\]

`A_k`, `D_k`, POP-H row errors, and cosine similarity are post-fit diagnostics and may not tune `mu`.

## 4. Remaining small-sample variance question

The primary asymptotic variance remains the unsmoothed cluster-score sandwich already specified in Section 8 of the base method specification. No scalar finite-sample multiplier is authorised.

For TRUE-SUPPORT diagnostics, decompose

\[
T_k=\sqrt n(\widetilde\beta_k-\beta_k^\star)
=L_k+R_k,
\]

where

\[
L_k=-n^{-1/2}\sum_i\omega_k^{\mathrm{pop}\,\mathsf T}g_i^\star.
\]

Report `Var(T_k)`, `Var(L_k)`, `Var(R_k)`, and `2Cov(L_k,R_k)`, together with a variance ladder comparing population first-order variance, target-score/sample covariance, fitted-score/sample covariance, and fitted-score/sample-exact-inverse covariance.

Also compute a delete-one-cluster jackknife and low-dimensional profile KC/CR2-style and MD/CR3-style leverage corrections as diagnostics. These are not promoted to the production high-dimensional variance until their mechanism and simulation behaviour are reviewed prospectively.

## 5. Pilot calibration regime

P01--P04 remain mandatory small-sample stress/diagnostic cells. Add:

```text
P05: n=200, p=500, s=5, tau=0.50, q=1, baseline Gaussian random intercept
P06: n=400, p=500, s=5, tau=0.50, q=1, same DGP
```

The existing `[0.80,1.20]` `mean(SE)/empirical SD` diagnostic and broad `[0.88,0.98]` coverage diagnostic are not relaxed. The hard first-order calibration assessment is made at P05; P06 is used to verify scaling. P01--P04 remain visible and retain correctness, identity, convergence, failure-accounting and mechanism-diagnostic requirements.

If P05 fails after the round-two precision and variance mechanism work, the final registry remains blocked and another theory review is required.

## 6. Numerical solver

HiGHS may replace CLARABEL as the production solver for the same Dantzig linear program only after a prospective parity suite on at least 20 representative rows verifies feasibility/infeasibility agreement, residual `<= mu+1e-6`, `l1` objective agreement within relative tolerance `1e-5` (plus a small absolute tolerance near zero), one-step output agreement within a committed tolerance, and deterministic repeatability. CLARABEL remains the reference audit solver.

## 7. Asset definitions

Profile-target assets continue to approximate the unsmoothed profile target using their own small numerical target-approximation bandwidth and may be reused under the approved dependency hash. POP-H assets use

\[
h_{\mathrm{analysis}}=c_h n_{\mathrm{analysis}}^{-3/10},
\]

with the large population sample controlling Monte Carlo accuracy only. The fields `n_analysis`, `h_analysis`, `sigma0_pop`, and `Sigma0_population` remain required audit fields.
