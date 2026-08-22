# Method specification — round-three pilot amendment

**Status:** authoritative amendment pending incorporation into `docs/METHOD_SPECIFICATION.md`  
**Date:** 2026-08-23  
**Source decision:** `results/preflight/PILOT_V2_THEORY_DECISION_ROUND3.md`  
**Scope:** first-stage penalty/loadings, first-stage numerical acceptance, lambda--mu sequencing, and pilot-gate interpretation

This amendment supersedes the raw rule `lambda_beta=sqrt(log p/n)` with unit penalty factors. It does **not** alter the exact profile score, corrected effective Hessian, target definition, one-step formula, primary unsmoothed meat, analysis bandwidth, or the round-two Dantzig program.

## 1. Cluster self-normalised first-stage penalty

Let `P` be the penalised coordinate set, `U` the always-included set and `p_P=|P|`. At a provisional coefficient vector `b`, define the exact smoothed cluster score contribution

\[
s_i(b)
=-m_i^{-1}X_i^{\mathsf T}
\psi_{\tau,h}\{r_i(b)\}.
\]

For `j in P`, define

\[
\bar s_j(b)=n^{-1}\sum_i s_{ij}(b),
\qquad
\widehat\ell_j(b)
=\left[n^{-1}\sum_i
\{s_{ij}(b)-\bar s_j(b)\}^2\right]^{1/2}.
\]

Freeze

\[
\alpha_{\lambda,n}=0.10/\log\{\max(n,3)\},
\]

\[
q_{\lambda,n}
=\Phi^{-1}\left(1-\frac{\alpha_{\lambda,n}}{2p_P}\right),
\qquad
\lambda_{0,n}=1.10\,q_{\lambda,n}/\sqrt n.
\]

Because the loading uses the actual `psi_{tau,h}` cluster score, do not insert an additional factor `sqrt{tau(1-tau)}`.

### 1.1 Exact two-pass algorithm

1. Set `b^(0)=0`; profile the nuisance effects and compute `ell_j^(0)=ell_j(b^(0))`.
2. Fit a preliminary estimator with coordinate penalty
   \[
   p_j^{(0)}=\lambda_{0,n}\widehat\ell_j^{(0)}
   \]
   on `P` and zero penalty on `U`. Call it `b^(1)`.
3. Compute `ell_j^(1)=ell_j(b^(1))` and
   \[
   \widehat\ell_j^{\mathrm{final}}
   =\max\{\widehat\ell_j^{(0)},\widehat\ell_j^{(1)}\}.
   \]
4. Warm-start at `b^(1)` and refit with
   \[
   \boxed{
   p_j=\lambda_{0,n}\widehat\ell_j^{\mathrm{final}}
   }
   \]
   on `P` and zero penalty on `U`.

Use exactly one loading update. No coefficient truth, support truth, empirical bias, coverage, p-value, or test-set outcome may enter the penalty.

### 1.2 Current optimiser mapping

Use

```text
lambda_beta = lambda_0_n
penalty_factor_j = ell_j_final * base_penalty_factor_j
```

where `base_penalty_factor_j=0` on `U` and `1` on ordinary penalised coordinates.

### 1.3 Degenerate loadings

Let `ell_ref` be the median positive loading among `P`. A retained coordinate with non-finite loading or

```text
ell_j < 1e-8 * ell_ref
```

causes `lambda_loading_degenerate`. Do not apply an arbitrary floor. Outcome-blind zero-variance design columns must be removed before this stage.

## 2. Lambda multipliers

The field `lambda_beta_multipliers` is reinterpreted as a prespecified sensitivity multiplier `a` applied to the calibrated coordinate penalties:

\[
p_j(a)=a\lambda_{0,n}\widehat\ell_j^{\mathrm{final}}.
\]

The primary method uses `a=1`. Alternative values occur only in separately registered sensitivity cells. The primary runner must not select `a` by coverage, truth, or an unrecorded search.

## 3. First-stage KKT and numerical acceptance

For final effective coordinate penalty `p_j`, define

\[
r_j=
\begin{cases}
|g_j+p_j\operatorname{sign}(\widehat\beta_j)|,
&j\in P,\ |\widehat\beta_j|>10^{-10},\\
\max(0,|g_j|-p_j),
&j\in P,\ |\widehat\beta_j|\le10^{-10},\\
|g_j|,&j\in U.
\end{cases}
\]

Let `p_ref` be the median positive `p_j` and define

\[
R_{\mathrm{KKT}}
=\max\left\{
\max_{j\in P}r_j/p_j,
\max_{j\in U}r_j/p_{\mathrm{ref}}
\right\}.
\]

A final first-stage fit is accepted only when

```text
R_KKT <= 1e-3
max_nuisance_gradient <= 1e-7
last beta_change <= 1e-7 * max(1, ||beta_hat||_inf)
all coefficients and objective values are finite
```

Reference solver controls:

```text
max_iter >= 2000
max_backtrack >= 50
beta_tol = 1e-7
```

Failure to meet the conditions is a `penalised_fit` failure. Objective stabilisation or a nonzero support is not a substitute for KKT acceptance.

## 4. First-stage theoretical rate

Let `beta_h_star` be the smoothed population profile minimiser. Under standardised score domination, uniformly comparable score loadings, profile restricted strong convexity, sparsity/moment conditions and numerical KKT error `o_p(lambda_0_n)`, retain

\[
\|\widehat\beta-\beta_h^\star\|_1
=O_p\left\{s\sqrt{\log p/n}\right\},
\]

\[
\|\widehat\beta-\beta_h^\star\|_2
=O_p\left\{\sqrt{s\log p/n}\right\}.
\]

The required score-domination event is

\[
\max_{j\in P}
\frac{|g_{h,j}(\beta_h^\star)|}{\ell_j^\star}
\leq\lambda_{0,n}/(1+\eta)
\]

for fixed `eta>0`. For the unsmoothed target,

\[
\|\widehat\beta-\beta^\star\|_1
=O_p\left\{s\sqrt{\log p/n}
+\|\beta_h^\star-\beta^\star\|_1\right\}.
\]

Under the manuscript target-approximation condition this becomes `O_p{s sqrt(log p/n)+s h^2}`.

The expression `s sqrt(log p/n)` is a stochastic rate, not a unit-constant finite-sample threshold. Simulations must nevertheless show shrinking error across the registered `n` sequence; an order-one empty fit is outside the intended regime.

## 5. Dantzig rule after lambda calibration

Round-three makes **no change** to the round-two Dantzig anchor, grid or non-truth selector:

\[
\mu(c)=c\left\{\sqrt{\log p/(nh)}+h^2\right\},
\qquad
c\in\{0.02,0.05,0.10,0.25,0.50,1,2,4\}.
\]

Use deterministic two-fold cluster CV for `n<200` and four-fold cluster CV for `n>=200`; minimise the mean held-out inverse defect and choose the smaller `mu` in a numerical tie.

The precision stage is conditional on an accepted final first stage. Recompute nuisance effects, per-cluster Hessian contributions and the full effective Hessian at the accepted final `beta_hat`; then select `mu`. A Hessian inherited from `beta=0`, the preliminary fit or a stale object is invalid.

The fact that `sqrt(n) mu` does not vanish is not the relevant criterion. The debiasing defect is controlled by

\[
\sqrt n\,\mu\|\widehat\beta-\beta^\star\|_1.
\]

Reopen the `mu` rate only if, after lambda repair, the actual `A_k`, POP-H row errors, total Bahadur remainder and P05--P06 scaling remain poor jointly. `D_k` alone is not a hard gate.

## 6. Oracle variance and pilot gate

The primary variance remains the corrected smoothed profile bread with the unsmoothed fitted cluster-score meat. No finite-sample scalar multiplier is authorised.

TRUE-SUPPORT/POP-H coverage is mechanism evidence, not an independent hard freeze gate. The observed `n=200` oracle coverage `0.87--0.92`, with negligible bias and `SE/SD=0.80--0.95`, is reported as mild finite-sample undercoverage and does not reopen the variance formula.

The practical proposed-method P05 gate retains

```text
SE/SD in [0.80,1.20]
coverage in [0.88,0.98]
```

and P06 remains the scaling confirmation. KC/CR2, MD/CR3 and delete-one-cluster jackknife remain diagnostics/sensitivities only.

## 7. Versioning and invalidation

The lambda change invalidates every practical P01--P06 fit for freeze purposes. Create a new versioned pilot registry and new run IDs; preserve old outputs. Profile-target and POP-H assets are unchanged mathematically and may be reused under a matching dependency hash because analysis lambda is not one of their dependencies. Practical fitted Hessians, selected precision rows, intervals and summaries must be regenerated.

## 8. Required audit fields

Every proposed-method fit stores at least:

```text
lambda_rule
lambda_alpha
lambda_normal_quantile
lambda_safety_constant
lambda_base
lambda_loading_pass0_min/median/max
lambda_loading_pass1_min/median/max
lambda_coordinate_min/median/max
zero_profile_score_max
zero_kkt_ratio_max
preliminary_kkt_normalized
final_kkt_normalized
final_kkt_absolute
first_stage_iterations
first_stage_beta_change
first_stage_nonzero_count
```

Simulation-only first-stage `l1/l2` errors are saved as diagnostics and never used for tuning.
