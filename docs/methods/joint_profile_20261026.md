# Joint profile estimation and inference

The current manuscript method is `JOINT-PROFILE`, implemented by the research estimator `joint_profile_design_20261026`. Its implementation guide is `research/score_validity_20261026/README_zh.md`. The population and selection arguments are stated in the manuscript and its Supporting Information.

## Observations and target

Independent clusters are `D_i=(Y_i,X_i,Z_i,m_i)`, for `i=1,...,n`. Cluster sizes satisfy `1 <= m_i <= m_max` in the asymptotic model, the fixed effect dimension may exceed `n`, and the nuisance dimension is fixed. The working nuisance matrix `Lambda` is positive definite and fixed during coordinate inference.

Define

\[
 \ell_{i,b}(\beta,\gamma)=m_i^{-1}\sum_j
 \rho_{\tau,b}(Y_{ij}-X_{ij}^T\beta-Z_{ij}^T\gamma)
 +\tfrac12\gamma^T\Lambda\gamma,
 \qquad \ell_{i,b}^{prof}(\beta)=\min_\gamma\ell_{i,b}(\beta,\gamma).
\]

The target is `beta_star=argmin E ell_i,0^prof(beta)`. It is defined at bandwidth zero. Structural slopes are valid simulation targets only under the corresponding population orthogonality condition; other settings use a separately computed profile target.

## Complete profile geometry

At the nuisance minimiser, write `r_i=Y_i-X_i beta-Z_i gamma_i,b(beta)`,
`W_i=diag(K_b(r_i))/m_i`, and `B_i=Z_i^T W_i Z_i+Lambda`. Then

\[
 g_{i,b}(\beta)=-X_i^T\psi_{\tau,b}(r_i)/m_i,
 \qquad H_{i,b}(\beta)=X_i^TW_iX_i-X_i^TW_iZ_iB_i^{-1}Z_i^TW_iX_i.
\]

The equivalent residualised-design expression contains `A_i^T Lambda A_i`. This term remains part of every curvature calculation. Each constituent Hessian is positive semidefinite. The reference kernel is Epanechnikov.

For a single observation, let `x=X_i^T`, `z=Z_i^T`,
`a=z^T Lambda^(-1) z` and `w=K_b(r_i)`.
Nuisance stationarity gives `gamma=Lambda^(-1) z psi_b(r_i)` and
`grad_beta r_i=-x/(1+a*w)`. Its profile gradient is `-x psi_b(r_i)` and
its Hessian is `x x^T*w/(1+a*w)`. The nuisance criterion remains strongly
convex because `Lambda` is positive definite. In the population expansion,
the active dual chart has either zero or one interior coordinate. The
existing section-density bounds and finite chart count apply to both;
the empirical envelopes and cluster variance calculations use the same
within-cluster average. Single observations therefore use the same target,
paired coefficients and influence formula. Supporting Information states
these identities, and the actual singleton cows are checked by
`research/milk_application_20261031/check_singleton_geometry.R`.

## Selection on independent clusters

For data dependent selection, use two outer folds of whole clusters. The complementary sample selects the model for an estimation fold. Three inner folds choose the outcome search size by exact heldout unsmoothed profile risk. Search steps use a weighted orthogonal basis and the complete profile dual score. Each selected joint convex fit uses prediction ridge `1/n_training` and a primal dual gap certificate.

For an inner sample with `r` clusters and `d_0` always included coordinates, the largest number of additions is

\[
 \min\{p-d_0,[\min\{\lfloor r/\log r\rfloor,
             \lceil(\log r)^3\rceil\}-d_0]_+\}.
\]

The selected supports from the three inner fits are combined. Target-specific profile graphs are learned from positive curvature factors at `h` and `h/2`, using geometric L1 penalty ratios and heldout quadratic loss. Their graph union uses capacity `B(r)=min(p,floor(r/log(r)),ceil(log(r)^3))`, before subtracting always-included coordinates. A separate nodewise regression uses the cluster weighted design without outcome weights. Its five cluster fold prediction loss selects a ratio from 24 values between one and `min(0.05,n_selection^(-1/2))`. Its retained support is bounded by the smallest inner `B(r)` and ordered by absolute standardised coefficients. The final model contains the outcome and profile graph supports, the design regression support and the requested target.

The theoretical target approximation is controlled by comparison coefficients and contrasts supported on this union. Extra design-selected variables preserve those comparison vectors. Target identities are fixed before observing the corresponding estimation outcomes.

## Paired coefficient estimation

On estimation fold `f`, with selected set `A_fk`, fit the complete criterion at both bandwidths `b=h,h/2` with the same selected set:

\[
 \widehat\beta_{f,b,A}\in\arg\min_{\beta_A}
 \left\{n_f^{-1}\sum_{i\in I_f}\ell_{i,b}^{prof}(\beta_A)
       +\frac{1}{2n_f^2}\beta_A^T\widehat\Sigma_{f,A}\beta_A\right\},
 \quad\widehat\Sigma_{f,A}=n_f^{-1}\sum_{i\in I_f}X_{i,A}^TX_{i,A}/m_i.
\]

The numerical prediction ridge is included in the objective and its derivative. The point estimate is

\[
 \widehat\theta_k=\sum_f(n_f/n)
            (2\widehat\beta_{f,h/2,A_{fk},k}-\widehat\beta_{f,h,A_{fk},k}).
\]

The default bandwidth is `h=n^(-3/10)`. Prespecified bandwidth sensitivity multipliers remain explicit. The complete profile population gradient can have a first order bandwidth term, and the paired coefficient construction leaves a target approximation of order `h^2` under the density and local curvature conditions.

If a design is fixed independently of its estimation outcomes, use one estimation sample with `A` fixed. This is the one-fold form of the same definition. It applies to the prespecified clinical and immune module model, and to a selected gene model evaluated on a separate group of subjects.

## Coordinate uncertainty

Let `omega_f,b` be the target row of the inverse selected Hessian, including its numerical ridge. Using the complete cluster score at the corresponding fitted coefficient, form

\[
 V_{ik}=-2\widehat\omega_{f,h/2}^Tg_{i,h/2}
           +\widehat\omega_{f,h}^Tg_{i,h},
 \qquad \widehat s_k^2=n^{-2}\sum_f\sum_{i\in I_f}(V_{ik}-\bar V_{fk})^2.
\]

The normal interval is `theta_hat +/- z_(1-alpha/2)*se`. Each cluster combines both bandwidth contributions before variance calculation. The numerical ridge has an asymptotically negligible influence contribution. For a prespecified finite family, Bonferroni allocation uses the declared family size.

Full workflow subject deletion repeats feature construction and fitting at the fixed original bandwidth. Its standard error is computed from the resulting coefficient vectors. Profile and full workflow intervals retain their distinct sampling definitions.

Exploratory gene designs use the frozen independent subject roles. Three forward additions to the clinical basis select the three genes, with `Lambda=I`. The separate estimation subjects fit the complete clinical and selected-gene design at two bandwidths. The interval target is that role's selected profile parameter. The role-specific intervals remain separate when the covariates or their measurement transformations differ.

## Population curvature

The curvature moment argument uses the active dual charts for any fixed nuisance dimension. On a chart with interior dual coordinates `A`, let `C_A=Z_A Lambda^{-1} Z_A^T/m`, `D_A=diag(k'(s_A))` and `Q_A=C_A+hD_A`. Its Jacobian is `det(Q_A)`, and its complete profile Hessian is `X_A^T Q_A^{-1} X_A/m`. If the kernel height is bounded by `K_max`, then

\[
 \det(Q_A)\|Q_A^{-1}\|_{op}^{r}
 \leq (K_{max}/h)^{r-1}\operatorname{tr}\{\operatorname{adj}(Q_A)\}.
\]

Each cofactor is a finite sum of bounded principal minors of `C_A` multiplied by products of `h k'(s_j)`. For `h<=1` these are dominated by an integrable product `C product_j(1+k'(s_j))`. The conditional density section bounds and the design contrast moments therefore give `E|H_{i,h,kl}(beta)|^r <= r! C^r h^(1-r)`. The bound is uniform over bounded nuisance designs with a fixed positive-definite `Lambda`, including repeated design rows and rank-deficient nuisance designs. The corresponding sparse directional bound supplies the learning and variance arguments.

The executable chart check covers nuisance dimensions one through four, repeated and close design rows, non-diagonal working penalties and mixed interior/endpoint dual patterns. A complete three-dimensional nuisance fit checks the same inference interface.

## Prediction for new subjects

Prediction uses fixed effects estimated without any test-subject nuisance adjustment. An outer training sample determines all probe choices and scaling. Its inner subject folds choose an unsmoothed joint profile search model by fixed-part pinball loss. The candidate working penalties are `Lambda=0.5I,I,2I`, as declared in the application configuration. The search budget counts additional covariates, `min(p-d_0,floor(r/log r))`, while all prespecified clinical and module covariates remain included. The selected penalty and model size are fitted on the complete outer training sample. Prediction at a new design row is `X_new^T beta_prediction`.

This prediction rule has its own validation loss and does not alter the fixed working matrix used for coordinate intervals. Clinical and module scores retain the declared public membership and fold-local transformations.

The random-intercept predictor uses the joint primal-dual solver. The random-slope predictor uses Clarabel on the equivalent sparse quadratic programme in the selected weighted QR basis. Nuisance RMS scaling transforms both `Z` and `Lambda`, preserving the complete criterion. Both solvers return a box-feasible dual score and are checked by the same recomputed normalised primal-dual gap.

## Numerical and output requirements

The nuisance equations, coefficient gradients and primal dual gaps must meet their declared tolerances. Numerical failure terminates the attempted unit and preserves its failure record. Results retain selection subjects, estimation subjects, selected variables and the covariance contributions at both bandwidths. No adjustment is selected using simulation truth or empirical coverage.

The complete profile geometry and conditional target definitions are shared with `reference_profile_pair_20261022.md`. That file defines the separately versioned paired-score reference estimator used in archived comparisons. Current manuscript estimation and uncertainty follow the paired coefficient equations above.
