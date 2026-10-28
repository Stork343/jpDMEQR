# Joint Profile Inference with Conditional Covariate Selection

The estimator `ds_fit_profile` has method identifier `joint_profile_design_20261026`. It targets the complete unsmoothed population profile parameter

\[
 \beta^\star(\Lambda)=\arg\min_\beta E\min_\gamma
 \left\{m_i^{-1}\sum_j\rho_\tau(Y_{ij}-X_{ij}^T\beta-Z_{ij}^T\gamma)
                    +\tfrac12\gamma^T\Lambda\gamma\right\}.
\]

`Lambda` is fixed and positive definite. Complete clusters are the independent observations. The implementation supports a random intercept and a random intercept with a time slope.

## Selection and estimation samples

Two folds partition whole clusters. Each fold selects variables using its complementary sample and fits the selected model using its own observations. Every cluster contributes to estimation once. Both estimation models explicitly include each requested target coordinate.

Within each selection sample, three inner folds construct dual guided orthogonal search paths. Every step scans all fixed effect columns, adds a variable and jointly optimises the selected coefficients and cluster nuisance effects. Weighted QR provides the orthogonal coordinates. A prediction ridge of `1/n_training` supplies a computable primal dual certificate. Exact heldout unsmoothed profile risk selects the model size, and the three selected inner supports are combined.

The inner search budget is the smaller of `floor(n_training/log(n_training))` and `ceiling(log(n_training)^3)`, after accounting for always included variables. The first budget binds at the evaluated sample sizes. Zero inclusion indicators specify always included covariates.

## Selecting covariates for each target

The outcome search is supplemented by the variables needed to adjust a target coefficient. Positive profile Schur factors at `h` and `h/2` define nodewise L1 paths. Their geometric penalty ratios decrease by `2^(-1/4)` from one to at most `n_training^(-1/2)`. At each ratio, the outcome support and both graph supports are combined, followed by a positive quadratic refit. Summed inner heldout quadratic losses select the ratio. The resulting supports are combined across the inner folds.

A separate nodewise regression uses only the cluster weighted design factor `X_ij/sqrt(m_i)`. Five cluster folds select its penalty by prediction error. Its ratio grid has 24 log spaced values from one to `min(0.05,n_selection^(-1/2))`. This regression identifies conditional covariate associations without estimating outcome density weights. Its selected variables are included in the estimation model. At most the inner search budget is retained, ordered by the absolute standardised regression coefficients. The graph support and its tuning losses are recorded.

The final selected set contains the outcome and profile supports, the design regression support, and the target coordinate. The design regression is particularly useful when a covariate contributes little to outcome prediction but is needed to adjust the target coefficient. Its inclusion preserves the population approximation supplied by the profile support.

## Joint profile fitting on the estimation sample

For target `k` and estimation fold `f`, let `A_fk` be its selected set. At each bandwidth `b` in `{h,h/2}`, jointly fit all its coefficients by minimising

\[
 \frac1{n_f}\sum_{i\in I_f}\ell_b(D_i,\beta_A)
       +\frac{1}{2n_f^2}\beta_A^T\widehat\Sigma_{f,A}\beta_A,
 \qquad \beta_{A^c}=0,
\]

where `ell_b` is the complete nuisance profile and
`Sigma_hat_f,A=n_f^(-1)sum_i X_i,A^T X_i,A/m_i`.

The small prediction ridge ensures positive numerical systems and vanishes faster than the inferential scale. The bandwidth is `h=n^(-3/10)`, using the total cluster count. Newton updates use the complete profile gradient and Schur curvature. Armijo backtracking checks the complete objective. Both the gradient and coefficient step must meet their numerical tolerances. An iteration limit raises an error.

Let `theta_f,h` and `theta_f,half` be the fitted target coefficients. The estimate is

\[
 \widehat\theta_k=\sum_f\frac{n_f}{n}
                    (2\widehat\theta_{f,h/2}-\widehat\theta_{f,h}).
\]

The complete profile population expansion can contain a term linear in bandwidth. The same selected model and prediction ridge are used in both fits, so the paired coefficient estimate removes that term and leaves an approximation of order `h^2`.

## Matched cluster uncertainty

At each fitted bandwidth, solve the selected positive Hessian system for the target row `omega_f,b`. Use the original fixed effect design in the complete cluster score `g_i,b=-X_i,A^T psi_i,b/m_i`. The contribution for uncertainty is

\[
 V_{ik}=-2\widehat\omega_{f,h/2}^Tg_{i,h/2}
                   +\widehat\omega_{f,h}^Tg_{i,h}.
\]

The covariance between bandwidth fits is retained by combining both contributions within the same cluster. With centring within each estimation fold,

\[
 \widehat s_k^2=n^{-2}\sum_f\sum_{i\in I_f}
                       (V_{ik}-\overline V_{fk})^2.
\]

The confidence interval is `theta_hat +/- qnorm((1+confidence)/2)*se`. The default confidence level is 0.95. The fixed effect numerical ridge contributes an asymptotically negligible term, as specified in `JOINT_THEORY.md`.

## Computational structure

All high dimensional scans occur during selection. Estimation then uses only the selected columns. The coefficient Hessian has the dimension of that selected model, while every cluster nuisance system has dimension `q`. The method avoids constructing or inverting a complete `p` by `p` matrix.

`small_profile.R` implements the selected convex solver and paired coefficient fit. `geometry.cpp` computes the Schur curvature in compiled code. The compiled nuisance engine is the package's existing complete profile solver. The unsmoothed outcome search uses the exact finite cluster solver and joint primal dual algorithm documented in `positive_curvature_20261025`.

Every result retains its training and estimation subjects, selected supports and tuning losses. Numerical gradients and the two bandwidth covariance contributions are stored. Simulation requests and failed attempts remain explicit records. A compilation cache lock allows independent R workers to start safely.

## Statistical guarantee

Outcome risk approximation and profile contrast approximation control the difference between the selected model target and the full profile target. Adding the design regression support preserves those bounds. Independent estimation then yields a paired profile influence expansion and a consistent cluster variance under the conditions in `JOINT_THEORY.md`.

For Gaussian covariates independent of the cluster errors and nuisance effects, retaining the conditional neighbours of the target coordinate gives an additional exact population result: the selected model target coefficient equals the full target even when other outcome coefficients are omitted. The corresponding proposition explains the design regression and is checked directly in the simulation target audit.
