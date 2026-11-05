# Conditional quantile estimation with repeated clinical measurements

## Scientific parameter

Patient i contributes observations j=1,...,m_i. Patients are independent and m_i is bounded. The design includes an intercept, a prespecified clinical vector and high dimensional measurements. The prediction target is Q_tau(Y_ij | W_ij)=W_ij^T b_tau. The parameter b_tau is independent of an optimisation penalty or a working covariance.

If the linear conditional quantile model is an approximation, CQ-ENET has the explicit population projection

\[
b_\tau^\dagger=\arg\min_b E\left[m_i^{-1}\sum_j\rho_\tau(Y_{ij}-W_{ij}^Tb)\right].
\]

This definition gives every patient equal weight. The conditional quantile interpretation requires the linear model. The correlation assisted score below uses the stronger condition E[u_i(b_tau) | W_i,m_i]=0. This condition permits cross visit instruments. In particular, subsequent outcomes used as regressors for other visits require a separate justification of this conditional moment.

## Estimation over all candidate measurements

Write b=(a,b_P), where a contains the prespecified clinical coefficients. For a training sample of n patients, CQ-ENET minimises

\[
L_n(b)+\lambda\{\alpha\|b_P\|_1+(1-\alpha)\|b_P\|_2^2/2\},
\quad L_n(b)=n^{-1}\sum_i m_i^{-1}\sum_j\rho_\tau(Y_{ij}-W_{ij}^Tb).
\]

Every clinical column is retained through a basis for its full estimable span. Exact linear dependence is removed with a training weighted QR decomposition, and the retained basis must reconstruct the requested clinical design. The protein coefficients may be sparse or distributed across many measurements. No bound of floor(n/log n) is imposed on the number of nonzero prediction coefficients.

All imputation and scaling parameters are learned on the corresponding training patients. For N training rows, multiply each row of W and Y by N/(n m_i). Positive homogeneity of the check loss makes the unweighted mean check loss of these transformed rows exactly L_n. hqreg_raw provides an elastic net path on this transformed design with intercept=FALSE. The design already contains the transformed intercept. This operation leaves the coefficient penalty unchanged.

The selected model is fitted to the declared nonsmooth objective with Clarabel. The implementation recomputes the objective, stationarity, residual constraints and dual box restrictions. hqreg uses a Huber approximation along its path; the path objective and exact selected objective are recorded separately. Numerical tolerances and full path completion are checked before results are accepted.

Inner patient folds choose alpha and the relative lambda index using patient weighted pinball loss. The selected model is fitted on all outer training patients. Predictions are evaluated at the retained patient fold. The same preprocessing and evaluation units are used for comparisons with archived prediction methods. Raw predictions and the prespecified nonnegative sorted quantiles are both retained.

## Orthogonal score for a coefficient

For a chosen coordinate, write W_i b=d_i theta+X_i beta. Let u_ij= tau-I(Y_ij<=W_ij b), F_i=diag(f_ij), where f_ij is the conditional density at the target quantile, and C_i is a positive definite working covariance of u_i. Define

\[
\gamma=\arg\min_g E[(d_i-X_i g)^T F_i C_i^{-1}F_i(d_i-X_i g)],
\quad r_i=d_i-X_i\gamma.
\]

The estimating score and sensitivity are

\[
\phi_i(\theta,\beta,\gamma,F,C)=r_i^TF_iC_i^{-1}u_i,
\qquad J=E[r_i^TF_iC_i^{-1}F_i d_i].
\]

The projection normal equation gives E[r_i^TF_i C_i^{-1}F_iX_i]=0 and J=E[r_i^TF_iC_i^{-1}F_i r_i]. Hence the derivative of E phi_i with respect to beta vanishes at the true quantile coefficients. The derivative with respect to gamma is zero because E[u_i | W_i,m_i]=0. The same conditional moment makes the derivatives with respect to F and C zero when the parameter and the other nuisances are at their population values.

Training patients determine the quantile fit, density estimates, covariance and projection. The complementary patient fold supplies scores. For fold f with training estimate theta_f, use

\[
\widehat\theta=\sum_f\frac{n_f}{n}\left\{\widehat\theta_f+
\frac1{n_f}\sum_{i\in I_f}\frac{\widehat\phi_i}{\widehat J_f}\right\},
\quad
\widehat s^2=n^{-2}\sum_f\sum_{i\in I_f}
\left(\widehat\phi_i/\widehat J_f-\overline{\widehat\phi/\widehat J}_f\right)^2.
\]

A fold estimate and its influence contribution are transformed to the original outcome and target variable units before aggregation. The complete quantile fit is used for u_i. The correction uses a single projection for the requested coefficient and never requires a full p by p inverse matrix.

The computational working covariance is exchangeable, C_i=tau(1-tau){(1-rho)I+rho 11^T}. Estimate rho from pair products of out of fold quantile scores in the training sample. Projection to [0,0.95] keeps all bounded cluster sizes positive definite. This is a specified working family; negative dependence remains admissible for consistency under the conditional moment, with rho=0 as the fitted working choice. Compare rho=0 and the estimated rho using the same folds and nuisance prediction fits.

Density estimation uses fitted quantiles at tau-delta and tau+delta. The reciprocal spacing estimate is 2 delta / (q_(tau+delta)-q_(tau-delta)), with a positive spacing bound determined in training units. For asymptotic consistency delta must tend to zero and the spacing estimation error must be smaller than delta. Fixed finite sample spacing is a working approximation whose effect is assessed numerically.

The projection uses the transformed design blocks C_i^(-1/2) F_i X_i and response C_i^(-1/2) F_i d_i. Penalised least squares selects its regularisation by patient folds. The projection includes the remaining clinical design without penalty. Its empirical balance residual and sensitivity are recorded.

## Population argument and conditions

The check loss identity gives

\[
E[\rho_\tau(Y-q-t)-\rho_\tau(Y-q)\mid W]
=\int_0^t\{F_{Y|W}(q+s)-\tau\}\,ds.
\]

When q is the conditional tau quantile and its density is positive in a neighbourhood, this is uniquely minimised at t=0. Averaging with positive patient weights preserves the minimiser. Under the conditional moment, E phi_i=0 for any measurable positive definite working C_i and finite design measurable F_i. Consistency of the point root therefore does not require the working covariance to be the true covariance.

For the orthogonal expansion, require a conditional density with a bounded derivative near the quantile, bounded working covariance eigenvalues and a sensitivity separated from zero. Let a_n be the L2 prediction error of the preliminary conditional quantile, b_n the error of the learned instrument, and c_n its conditional derivative balance error. Sufficient remainder conditions are sqrt(n) a_n^2 -> 0, sqrt(n) a_n b_n -> 0 and sqrt(n) a_n c_n -> 0, together with consistent sensitivity and influence variance and a patient level Lindeberg condition. The instrument norm and moments must be controlled in these statements. Nuisance fitting uses a fixed number of patient folds.

Conditional on each training sample, the heldout empirical error of the estimated influence function is o_p(n^(-1/2)) when its L2 difference tends to zero. A Taylor expansion of the conditional distribution gives the quadratic prediction remainder. The projection identity removes the first order nuisance coefficient term; the estimated instrument and balance errors give the two product remainders. The independent patient central limit theorem then yields

\[
\sqrt n(\widehat\theta-\theta)=n^{-1/2}\sum_i J^{-1}\phi_i+o_p(1).
\]

These rate conditions are an explicit theorem interface. Deriving them from a growing dimensional elastic net model, learned quantile spacings and the chosen projection penalty requires additional assumptions and proof. The computational implementation does not by itself verify those rates in a clinical cohort.

When C_i equals the conditional score covariance, the displayed influence variance is J^(-1). In the linear instrument class generated by the covariates, it is the covariance weighted orthogonal projection for theta. A claim of semiparametric efficiency over all design measurable instruments would additionally require the corresponding projection over that larger nuisance tangent space.

## Research contribution to test

Elastic net quantile estimation and correlation weighted longitudinal estimating equations are established components. The research question is whether a patient cross fitted, target specific orthogonal construction can support interpretable quantile contrasts when prediction signals extend beyond a small selected subset, at a computational cost compatible with high dimensional clinical panels. The comparison must identify gains attributable to the projection and working correlation, in addition to gains from the more flexible preliminary estimator.

The closest work map is recorded separately with verified sources. Newness and improved performance remain research questions until those comparisons and the corresponding rate arguments are complete.
