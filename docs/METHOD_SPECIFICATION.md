# Authoritative method specification

## 1. Scope and asymptotic unit

We observe independent clusters

\[
D_i=\{(Y_{ij},X_{ij},Z_{ij}):j=1,\ldots,m_i\},\qquad i=1,\ldots,n,
\]

where `n` is the number of independent clusters, `m_i` is uniformly bounded in the primary asymptotic regime, `X_ij in R^p` is the fixed-effect design, and `Z_ij in R^q` is a low-dimensional cluster-effect design. The fixed-effect vector may be sparse with `p >> n`; `q` is fixed.

All empirical averages, stochastic rates and standard errors are indexed by `n`, not by the total row count treated as if all observations were independent. Within-cluster loss averaging is used so clusters do not receive accidental weight solely because they contain more observed visits.

## 2. Unsmoothed regularised profile target

For quantile level `tau in (0,1)`, define the check loss

\[
\rho_\tau(u)=u\{\tau-I(u<0)\}.
\]

Let `Lambda` be a fixed positive-definite `q x q` working penalty matrix. For a realised cluster define

\[
q_i(\beta,\gamma;D_i)
=\frac{1}{m_i}\sum_{j=1}^{m_i}
\rho_\tau(Y_{ij}-X_{ij}^{\mathsf T}\beta-Z_{ij}^{\mathsf T}\gamma)
+\frac{1}{2}\gamma^{\mathsf T}\Lambda\gamma.
\]

The cluster-specific profile map is

\[
\gamma_i^\dagger(\beta;D_i)
\in\arg\min_{\gamma\in\mathbb R^q}q_i(\beta,\gamma;D_i),
\]

and the unsmoothed population profile risk is

\[
Q(\beta)=E\left[q_i\{\beta,\gamma_i^\dagger(\beta;D_i);D_i\}\right].
\]

The primary inferential target is

\[
\beta^\star\in\arg\min_{\beta\in\mathbb R^p}Q(\beta).
\]

This is a regularised profile parameter. It depends on the working nuisance penalty `Lambda`, but it does not depend on the smoothing bandwidth. It need not equal a structural conditional-quantile coefficient under random-slope misspecification, informative cluster size, endogenous regressors, or other departures from the ideal location-shift model.

### 2.1 Structural coefficient versus profile target

When a simulation reports error or coverage relative to `beta0`, the scenario definition must justify `beta_star = beta0` for the coordinates being evaluated. A sufficient baseline design used in the simulations is:

- centred fixed-effect covariates independent of cluster effects and the random-effect design;
- no penalised fixed intercept;
- a correctly centred error distribution with conditional `tau`-quantile zero;
- a working nuisance design containing the generated random-effect directions;
- cluster size independent of covariates and outcomes.

Under misspecification, use a high-accuracy population Monte Carlo approximation to `beta_star` and propagate its numerical error into the coverage audit. Never silently substitute `beta0`.

## 3. Convolution smoothing

Let `K` be a symmetric probability density with zero first moment, finite second moment and bounded derivative. With `K_h(u)=K(u/h)/h`, define

\[
\rho_{\tau,h}(u)=(\rho_\tau*K_h)(u),\qquad
\psi_{\tau,h}(u)=\frac{\partial}{\partial u}\rho_{\tau,h}(u),\qquad
\phi_{\tau,h}(u)=\frac{\partial^2}{\partial u^2}\rho_{\tau,h}(u).
\]

The reference implementation uses the Epanechnikov kernel. Smoothing serves two roles: it makes the optimisation differentiable and provides an empirical local-curvature estimator. The paper targets `beta_star`, not the minimiser of a fixed-bandwidth smoothed risk.

For each cluster define

\[
q_{i,h}(\beta,\gamma)
=\frac{1}{m_i}\sum_{j=1}^{m_i}
\rho_{\tau,h}(Y_{ij}-X_{ij}^{\mathsf T}\beta-Z_{ij}^{\mathsf T}\gamma)
+\frac{1}{2}\gamma^{\mathsf T}\Lambda\gamma.
\]

The empirical profile objective is

\[
\widehat Q_h(\beta)
=\frac{1}{n}\sum_{i=1}^n
q_{i,h}\{\beta,\widehat\gamma_i(\beta)\},
\qquad
\widehat\gamma_i(\beta)=\arg\min_\gamma q_{i,h}(\beta,\gamma).
\]

The penalised estimator is

\[
\widehat\beta
\in\arg\min_\beta
\left\{
\widehat Q_h(\beta)
+\lambda_\beta\sum_{k=1}^p w_k|\beta_k|
\right\},
\]

where `w_k=0` for always-included covariates and `w_k=1` for penalised high-dimensional predictors unless adaptive weights are explicitly prespecified.

## 4. Exact profile calculus

For cluster `i`, write matrices `X_i` (`m_i x p`) and `Z_i` (`m_i x q`) and residual vector

\[
r_i(\beta)=Y_i-X_i\beta-Z_i\widehat\gamma_i(\beta).
\]

Define

\[
\Psi_i(\beta)=\psi_{\tau,h}\{r_i(\beta)\},
\qquad
W_i(\beta)=\frac{1}{m_i}\operatorname{diag}
\left[\phi_{\tau,h}\{r_{ij}(\beta)\}\right].
\]

### 4.1 Profile score

Because the nuisance minimiser solves the stationarity equation for the **complete** ridge-regularised cluster criterion, the envelope theorem gives

\[
\widehat g_h(\beta)
=\nabla\widehat Q_h(\beta)
=-\frac{1}{n}\sum_{i=1}^n
\frac{1}{m_i}X_i^{\mathsf T}\Psi_i(\beta).
\]

This score uses the original fixed-effect design `X_i`.

**Forbidden legacy formula:**

\[
-n^{-1}\sum_i m_i^{-1}\widetilde X_i^{\mathsf T}\Psi_i
\]

is not the profile score of the complete ridge criterion. Do not use it in the KKT expansion, one-step update, influence function, or sandwich variance.

### 4.2 Nuisance derivative

Define

\[
B_i(\beta)=Z_i^{\mathsf T}W_i(\beta)Z_i+\Lambda,
\]

\[
A_i(\beta)=B_i(\beta)^{-1}Z_i^{\mathsf T}W_i(\beta)X_i.
\]

Differentiating the nuisance stationarity equation gives

\[
\frac{\partial\widehat\gamma_i(\beta)}{\partial\beta^{\mathsf T}}
=-A_i(\beta).
\]

The sign convention follows residual `Y-X beta-Z gamma`; verify it numerically after any code change.

### 4.3 Effective profile Hessian

The exact Hessian is

\[
\widehat H_{\mathrm{eff}}(\beta)
=\frac{1}{n}\sum_{i=1}^n
\left[
X_i^{\mathsf T}W_iX_i
-X_i^{\mathsf T}W_iZ_i
B_i^{-1}
Z_i^{\mathsf T}W_iX_i
\right].
\]

Equivalently, with

\[
\widetilde X_i=X_i-Z_iA_i,
\]

we have the exact identity

\[
\widehat H_{\mathrm{eff}}(\beta)
=\frac{1}{n}\sum_{i=1}^n
\left[
\widetilde X_i^{\mathsf T}W_i\widetilde X_i
+A_i^{\mathsf T}\Lambda A_i
\right].
\]

The second term is positive semidefinite and is required. A unit test must confirm the maximum absolute difference between the two expressions is below `1e-8` in small deterministic examples and below a documented numerical tolerance in larger examples.

## 5. Penalised optimisation

The reference implementation uses a correctness-first profiled proximal-gradient algorithm:

1. Given `beta`, solve every low-dimensional nuisance subproblem to a strict gradient tolerance.
2. Compute the exact profile score.
3. Use the largest eigenvalue of the effective Hessian as an initial local Lipschitz estimate.
4. Apply a penalty-factor-aware proximal step.
5. Backtrack on the exact profiled composite objective.
6. Stop only when both parameter change and KKT residual meet tolerance.

A faster joint BCD implementation is permitted only after it is shown to return the same profiled solution and the same score/Hessian objects within tolerance.

### 5.1 KKT residual

For unpenalised coordinate `k`, the residual is `|g_k|`. For an active penalised coordinate it is

\[
|g_k+\lambda_\beta w_k\operatorname{sign}(\widehat\beta_k)|,
\]

and for an inactive coordinate it is

\[
\max\{0,|g_k|-\lambda_\beta w_k\}.
\]

The maximum coordinate residual must be recorded for every fit.

## 6. Row-wise precision estimation

For a prespecified target coordinate `k`, let `e_k` be the canonical vector. Estimate only the needed row/direction:

\[
\widehat\omega_k
\in\arg\min_{\omega\in\mathbb R^p}\|\omega\|_1
\quad\text{subject to}\quad
\|\widehat H_{\mathrm{eff}}\omega-e_k\|_\infty\leq\mu.
\]

The implementation must report:

- requested `mu`;
- achieved infinity-norm residual;
- solver and solver status;
- `||omega||_1` and `||omega||_2`;
- whether `mu` was enlarged for numerical feasibility.

If a feasibility grid is used, it must be deterministic and independent of the true coefficient and empirical coverage. The smallest feasible value in the frozen grid is used.

Full-matrix CLIME may be used for diagnostics, but coordinate-wise inference should not require estimating rows that are never used.

## 7. Debiased coordinate and influence representation

The one-step coordinate estimator is

\[
\widetilde\beta_k
=\widehat\beta_k
-\widehat\omega_k^{\mathsf T}
\widehat g_h(\widehat\beta).
\]

The intended expansion is

\[
\sqrt n(\widetilde\beta_k-\beta_k^\star)
=\frac{1}{\sqrt n}\sum_{i=1}^n
\xi_{ik}+o_p(1),
\]

where a sign-equivalent influence contribution is

\[
\xi_{ik}
=-\omega_k^{\star\mathsf T}g_i^\star,
\qquad
g_i^\star
=-m_i^{-1}X_i^{\mathsf T}\psi_\tau(r_i^\star).
\]

The sign is irrelevant for variance but must be consistent in any reported Bahadur remainder.

## 8. Cluster sandwich variance

### 8.1 Primary meat: unsmoothed quantile score

The inferential target is the **unsmoothed** regularised profile parameter, so the
primary Wald meat is built from the ordinary quantile score at the fitted profile
residuals (authoritative: `docs/PILOT_GATE_THEORY_DECISIONS.md`; the previous
smoothed-meat primary variance is superseded):

\[
\widehat g_i^{(0)}
=-\frac{1}{m_i}X_i^{\mathsf T}
\psi_\tau(\widehat r_i),
\qquad
\psi_\tau(u)=\tau-I(u<0),
\qquad
\overline g^{(0)}=n^{-1}\sum_i\widehat g_i^{(0)}.
\]

For coordinate `k`, compute

\[
\widehat\sigma_{0,k}^2
=\frac{1}{n}\sum_{i=1}^n
\left[\widehat\omega_k^{\mathsf T}
(\widehat g_i^{(0)}-\overline g^{(0)})\right]^2.
\]

The estimated finite-sample standard error is

\[
\widehat{\operatorname{se}}(\widetilde\beta_k)
=\widehat\sigma_{0,k}/\sqrt n.
\]

A nominal `(1-alpha)` Wald interval is

\[
\widetilde\beta_k
\pm z_{1-\alpha/2}\widehat\sigma_{0,k}/\sqrt n.
\]

The **bread remains the corrected smoothed effective Hessian** (via the
precision row). No scalar smoothing-deficit correction and no empirical SE
multiplier is permitted; the finite-`h` smoothing deficit is not a universal
additive correction in the clustered/profile model.

### 8.2 Diagnostic meat: smoothed score

The smoothed-score sandwich

\[
\widehat g_i^{(h)}
=-\frac{1}{m_i}X_i^{\mathsf T}\psi_{\tau,h}(\widehat r_i),
\qquad
\widehat\sigma_{h,k}^2
=\frac{1}{n}\sum_{i=1}^n
\left[\widehat\omega_k^{\mathsf T}
(\widehat g_i^{(h)}-\overline g^{(h)})\right]^2
\]

is retained **only as a diagnostic/sensitivity** (`se_smoothed`). Both meats
are reported side by side in the pilot and final outputs.

A degrees-of-freedom multiplier `n/(n-1)` may be reported as a prespecified sensitivity analysis; it cannot be chosen by coverage.

## 9. Tuning rules

### 9.1 Smoothing bandwidth

Primary inferential grid:

\[
h=c_h n^{-3/10},\qquad c_h\in\{0.75,1,1.25\}.
\]

The main specification uses `c_h=1`. This exponent satisfies both
`sqrt(n) h^2 -> 0` and `n h^3 -> infinity`; the former `n^{-1/3}` rule sits on
the boundary `n h^3 = 1` and is no longer admissible as the primary inferential
bandwidth. The full grid is used in the bandwidth module. Bandwidth selection
by validation loss is allowed for prediction analyses, but inferential
simulations report each prespecified bandwidth rather than choosing one by
coverage. The population-target approximation may continue to use its own
separately named small target-approximation bandwidth, which must be
independent of the analysis Dantzig tolerance.

### 9.2 Fixed-effect penalty

Use

\[
\lambda_\beta=c_\lambda\sqrt{\log(p)/n}
\]

with a frozen multiplier grid, for example `c_lambda in {0.25,0.5,1,2}`. The multiplier may be selected by cluster-level cross-validation on training clusters. The grid and tie rule must be frozen before simulations.

### 9.3 Nuisance penalty

The main setting is `Lambda=I_q`. Sensitivity values are `0.25 I_q`, `0.5 I_q`, `2 I_q` and `4 I_q`. In the paper, inference is conditional on the chosen working penalty. Estimating `Lambda` from the same data requires theory and code not yet included in the primary result.

### 9.4 Dantzig tolerance

A theory-motivated starting value is

\[
\mu=c_\mu\left\{\sqrt{\log(p)/(nh)}+h^2\right\},
\qquad h=c_h n^{-3/10}.
\]

The frozen candidate multiplier grid is

\[
c_\mu\in\{0.10,0.25,0.50,1,2,4\}.
\]

Feasibility of the constraint alone is not an adequate precision-direction
gate. The constant is selected **without simulation truth, bias or coverage**
by the frozen cluster-level inverse-Hessian cross-validation rule: partition
the independent clusters into four deterministic folds; solve the Dantzig row
on each training-fold Hessian at every candidate; score it on the held-out
fold with the inverse quadratic loss
`L(mu) = 0.5 omega' H_val omega - e_k' omega`; take the largest feasible
candidate within one standard error of the minimum mean loss; re-estimate the
row once on the full Hessian. Row-accuracy diagnostics (relative l1/l2 error
to the POP-H direction, cosine similarity, and the inverse-defect bound
`D_k = sqrt(n) delta_k ||beta_hat - beta_star||_1 / sigma_{0,k}^{pop}`) are
recorded as gates, never as tuning criteria.

## 10. Theory diagnostics required in simulations

For small/moderate dimensions where an oracle approximation can be computed, record:

1. target approximation `||beta_h_star-beta_star||`;
2. penalised estimation error;
3. `||H_hat-H_star||_max`;
4. Dantzig row error and constraint residual;
5. standard-error ratio `mean(SE)/empirical SD`;
6. scaled Bahadur remainder;
7. cluster-score normal-QQ diagnostics;
8. log-log slopes against `n` for principal error terms.

These quantities link the finite-sample experiment to the proof architecture. Prediction loss alone is insufficient.

## 11. Numerical identity tests

The geometry validation script must test all of the following over multiple seeds:

- `tau` in `{0.25,0.5,0.75}`;
- `q` in `{1,2}`;
- balanced and unbalanced cluster sizes;
- `h` multipliers `{0.75,1,1.25}`;
- non-diagonal positive-definite `Lambda`;
- finite-difference steps around `1e-5` with a sensitivity check.

Default acceptance thresholds:

- profile gradient maximum absolute error `< 2e-4`;
- profile Hessian maximum absolute error `< 3e-3`;
- Schur/identity maximum absolute error `< 1e-8`;
- nuisance gradient maximum absolute error `< 1e-7`;
- Dantzig residual `<= mu + 1e-6`.

Any threshold relaxation must be justified in a committed note and cannot be used to conceal a systematic sign or scaling error.

## 12. Interpretation boundary

The method provides coordinate-wise inference for the chosen regularised profile target. It does not automatically provide:

- causal effects;
- inference for an unregularised random-effect distribution;
- simultaneous inference over thousands of post-selected genes;
- FDR control after same-sample screening;
- validity under informative cluster size or cross-cluster dependence;
- validity when `Lambda` is estimated without accounting for its uncertainty.

These are explicit limitations, not implementation details to be hidden.
