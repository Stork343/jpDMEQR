# Cluster Profile L1 Reference Estimator

Method identifier: `crossfit_positive_curvature_20261025`.

This research specification implements the requested estimator redesign. The inferential target is the minimiser of the unsmoothed complete ridge profile population risk with fixed positive definite `Lambda`. Independent units are clusters. The numerical bandwidth does not define the target.

## Initial estimation

Let `D` contain the cluster weighted RMS of each fixed effect column, let `theta = D beta`, and use `X_s = X D^{-1}`. At each loading pass solve the joint convex problem

\[
 \min_{\theta,\gamma_1,\ldots,\gamma_n}
 \frac1n\sum_i\left\{\frac1{m_i}\sum_j
 \rho_\tau(Y_{ij}-X_{s,ij}^{T}\theta-Z_{ij}^{T}\gamma_i)
 +\frac12\gamma_i^T\Lambda\gamma_i\right\}
 +\sum_k a_k|\theta_k|+\frac{\|\theta\|_2^2}{2n}.
\]

The deterministic `1/n` initial ridge coefficient makes the complete finite sample optimisation strongly convex. It is recorded explicitly. It contributes order `1/n` to the initial population stationarity equation for bounded coefficients. The inferential score excludes this initial ridge and all L1 penalties.

Coordinate penalties are `a_k=lambda_base*loading_k*penalty_factor_k`. Zero penalty factors are supported. The reference fit uses always included fixed effects and the exact unsmoothed nuisance profile. Its centred full design cluster score standard deviations define fixed loadings. Set `lambda_max=max(abs(mean_score_j)/(loading_j*penalty_factor_j))` over penalised coordinates. Eight geometric ratios from one to `0.05` define the path.

Five deterministic folds contain whole clusters. Each training fold recomputes scaling, reference coefficients, loadings and `lambda_max`, and uses warm starts along its path. Validation sums exact unsmoothed heldout cluster profile losses with equal cluster weight. This criterion estimates the specified profile risk and uses each heldout outcome to optimise that cluster's nuisance adjustment. It is distinct from predicting a new patient's outcome. The ratio with smallest total validation loss is used for the complete sample coefficient fit. Candidate gaps use tolerance `1e-7`; the complete sample fit uses `1e-8`. Generated coefficients never enter selection.

The selected heldout fits also supply the inferential scores. For cluster `i` in fold `f`, evaluate the score and curvature at `beta_f`, fitted without that fold. Set `beta_bar=sum_f (n_f/n)*beta_f`. All these coefficient fits already exist in the validation paths. Statistical control must be uniform over the finite declared family of candidate procedures. The initial ridge coefficient is `1/n_T` when the training fold has `n_T` clusters.

Writing the objective multiplied by `n`, the dual variable obeys
`(tau-1)/m_i <= u_ij <= tau/m_i`. With `v = X_s^T u` and `z_i = Z_i^T u_i`, a feasible dual value is

\[
 d(u)=Y^Tu-\frac12\sum_k S(v_k,n a_k)^2
             -\frac12\sum_i z_i^T\Lambda^{-1}z_i,
\]

where `S` is soft thresholding. Subtracting this value from the complete primal objective gives a nonnegative duality gap. All dual iterates lie in the box. The gap is valid with zero L1 penalty factors and gives an actual optimisation certificate.

The implementation uses Algorithm 1 of Malitsky and Pock with linesearch. Matrix products use the original design and cluster blocks. No full fixed effect Hessian is needed. The two proximal operations are weighted soft thresholding divided by `1+step`, and `solve(I+step*Lambda, gamma+step*Z^T u)`. Dual updates are box projections. The linesearch tests the actual operator increment, so it does not require a global spectral decomposition.

## Score and positive curvature

At each cluster's heldout coefficient vector solve the nuisance profiles at `h=n^(-3/10)` and `h/2`, using the complete sample cluster count for this fixed bandwidth. Set

\[
 G_i^R=2g_{i,h/2}-g_{i,h},\qquad
 H_{i,h}=X_{s,i}^{T}L_{i,h}X_{s,i},
\]
\[
 L_{i,h}=W_i-W_iZ_i(Z_i^TW_iZ_i+\Lambda)^{-1}Z_i^TW_i,
 \qquad W_i=\operatorname{diag}\{K_h(r_{ij})\}/m_i.
\]

The positive semidefinite matrix `L_i,h` has dimension `m_i`. Its symmetric square root produces `F_i=L_i,h^(1/2) X_s,i`. Storing the combined rectangular factor `F` costs order `N*p`; `F^T F/n` equals the positive single bandwidth Hessian. The paired score is used in the correction and variance. Precision directions use the positive single bandwidth factor.

## Nodewise directions

For target `k`, standardise the columns of `F` by `d_j=sqrt(sum(F_j^2)/n)`. Regress column `k` on all other columns with a weighted L1 criterion. The objective is

\[
 \frac{1}{2n}\sum_i\|F_{i,k}/d_k-F_{i,-k}\operatorname{diag}(d_{-k}^{-1})a\|_2^2
 +\sum_{j\ne k}\lambda_{kj}|a_j|.
\]

The coordinate penalty is common within a nodewise regression. Twenty four geometric ratios from one to `0.05` define a path relative to the maximum absolute empirical predictor response product. Five cluster folds select the ratio using total heldout residual sum of squares. Each training fold defines its own anchor and uses its cluster count in the objective. These losses assess the conditional approximation of the cross fitted curvature factor. They are not prediction test errors. The full data direction uses the selected ratio and its full data anchor.

For residual `r_k`, set `t_k^2=sum(r_k^2)/n+sum(lambda_k*abs(a_k))`. Embed `c_k=(1,-a_k)` in the full coordinate order and set `omega_k=diag(d^{-1}) c_k/(d_k*t_k^2)`. The fitted nodewise KKT equations imply `(H_h omega_k)_k=1` and bound the other coordinates by their weighted penalties divided by `d_k*t_k^2`. The implementation checks both statements in the original standardised coordinates.

The rectangular design is passed to the published `glmnet` coordinate descent implementation. Its objective normalisation and penalty factor renormalisation are compensated explicitly. Full covariance matrices are not formed.

## Coordinate correction

\[
 \widetilde\theta_k=\overline\theta_k-\widehat\omega_k^T\overline G^R,
 \qquad
 \widehat{\mathrm{se}}(\widetilde\beta_k)^2
 =\frac{1}{n^2D_{kk}^2}\sum_i
 [\widehat\omega_k^T(G_i^R-\overline G^R)]^2,
 \qquad \widetilde\beta_k=\widetilde\theta_k/D_{kk}.
\]

All original fixed effect columns participate. The procedure does not construct a selected support refit. The requested coordinates need not have nonzero initial estimates. The two bandwidth scores are combined within each cluster before computing variance.

## Statistical derivation

Write `Delta_f=theta_f-theta_star`, `Delta_bar=sum_f(n_f/n)*Delta_f`, and let `H_hat` denote the positive curvature evaluated at the respective heldout coefficients. The exact identity is

\[
 \widetilde\theta_k-\theta_k^\star
 =-\widehat\omega_k^TG_n^R(\theta^\star)
 +(e_k-\widehat H\widehat\omega_k)^T\overline\Delta
 +\widehat\omega_k^T\left\{\widehat H\overline\Delta-
 \frac1n\sum_f\sum_{i\in I_f}[G_i^R(\theta_f)-G_i^R(\theta^\star)]\right\}.
\]

The final term includes the deliberate difference between positive curvature and the paired score derivative. Under the existing uniform profile expansion and curvature conditions it is bounded by
`O_p(||omega_hat||_1 * {(h+c_n)*e_n+e_n^2})`, where `e_n` bounds the initial L1 error and `c_n` bounds the uniform empirical Hessian error. Nodewise optimisation residuals contribute to the second term through the checked inverse defect `mu_n`.

A conditional asymptotic normality result follows if the population paired score is order `h^2`, the population inverse row is L1 bounded, its estimate has L1 error `r_n`, and

\[
 \sqrt n[\mu_n e_n+(h+c_n)e_n+e_n^2+h^2]\to0,
 \qquad r_n\sqrt{\log(2p)}\to0,
 \qquad e_n^2/h\to0.
\]

The limiting influence uses the inverse of the unsmoothed population curvature. The difference between that inverse and the single bandwidth population inverse is order `h` under bounded inverse norms. A cluster Lindeberg condition and convergence of the projected score variance complete the central limit theorem. Uniform score substitution and direction consistency give sandwich consistency.

The initial estimation result requires restricted Bregman curvature and score domination, including the initial ridge contribution and numerical error. It yields the candidate rate `e_n=O_p(s*sqrt(log(p)/n))` when the selected penalty satisfies these conditions and the numerical contribution is negligible. Primitive guarantees for the actual profile risk validation and nodewise validation rules require further analysis. The exact decomposition and a conditional inference theorem are developed in `THEORY.md`. Numerical coverage never determines a penalty, bandwidth or interval multiplier.

## Experimental definition

Numerical validation compares the initial solver with an independently formulated conic quadratic programme. Factor identities and nodewise KKT conditions are checked for random intercepts and intercept plus slope designs, including non diagonal `Lambda`. Scaling measurements use registered A02, A05 and A08 data generation with the seed retained in the output. Statistical assessment includes three quantiles and both nuisance dimensions. It reports each requested coordinate, convergence, variable selection, standard error ratios, coverage and Monte Carlo uncertainty.

Only observations and declared tuning rules enter the estimator. Generated coefficients are used by the evaluation code after estimation. Existing final registries and application results retain their own method identifiers.

## Sources and contribution

Malitsky and Pock provide the primal dual linesearch algorithm. Van de Geer and colleagues provide nodewise inverse construction. These components are attributed to their authors. The research contribution under investigation is the combination of exact unsmoothed joint fitting, positive profile curvature and a paired score correction, together with the curvature mismatch remainder and computable optimisation errors in clustered quantile inference. The literature search does not establish a priority claim.

- Malitsky and Pock 2018 [A First Order Primal Dual Algorithm with Linesearch](https://arxiv.org/abs/1608.08883)
- Van de Geer and colleagues 2014 [On Asymptotically Optimal Confidence Regions and Tests for High Dimensional Models](https://arxiv.org/abs/1303.0518)
- [glmnet documentation](https://glmnet.stanford.edu/articles/glmnet.html)
