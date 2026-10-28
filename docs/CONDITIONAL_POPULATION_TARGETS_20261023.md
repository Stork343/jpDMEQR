# Conditional population targets under covariate dependent scale

The profile parameter is the minimiser of the complete unsmoothed ridge-profile population risk. In B16 and B19, the registered scale is `exp(0.35*clip(X1,-3,3))`. Conditional quantile centring alone leaves its profile score mean unrestricted. These cases use a population target built by conditional Gaussian integration.

## Exact population reduction

The design has independent Gaussian AR1 rows. Its covariance is `Sigma_kl=rho^|k-l|`. Put `a_k=rho^(k-1)` and `C=Sigma-aa'`. Conditional on every first-coordinate value within a cluster,

\[
X_j=aX_{j1}+\eta_j,\qquad E\eta_j=0,\quad
\operatorname{Cov}(\eta_j)=C.
\]

The residual rows eta are independent across visits and independent of cluster size, nuisance design, cluster effects and error innovations. The generating fixed-effect signal is retained in full. At `beta=beta0+delta e1`, subtracting `X beta0` leaves the nuisance contribution, the X1-dependent error and `-delta X1`. The nuisance optimizer consequently depends on none of the eta coordinates. For either an ordinary envelope score or its paired version,

\[
E\{g_i(\beta^0+\delta e_1)\mid X_{i1},Z_i,m_i,b_i,\epsilon_i\}
=a\{-X_{i1}^{\mathsf T}s_i(\delta)/m_i\}.
\]

The scalar root therefore solves every population score coordinate. The within-cluster denominator stays inside the expectation, so random independent cluster sizes preserve this identity. The unsmoothed risk is convex. Positive conditional outcome densities and the ridge penalty give positive curvature at the scalar root and identify its full-dimensional minimum.

The numerical integration generates the marginal design containing all nonzero generating coordinates and three additional coordinates, up to the registered dimension. For B16 and B19 this is the first eleven coordinates of the full 1000-coordinate law. All remaining coordinates are integrated out using the displayed Gaussian conditional distribution. The object records this integration definition and sets `finite_full_design_empirical_root=FALSE`.

## Curvature and covariance

With `W_a=diag(K_a(r))/m`, write

\[
V_a=W_a-W_aZ(Z^{\mathsf T}W_aZ+\Lambda)^{-1}Z^{\mathsf T}W_a.
\]

For `V_R=2V_half-V_full`, conditional integration gives

\[
A_h=\alpha_h aa^{\mathsf T}+\beta_h C,\qquad
\alpha_h=E(X_1^{\mathsf T}V_RX_1),\quad
\beta_h=E\operatorname{tr}(V_R).
\]

Both coefficients are evaluated together with the residualised-design identity containing the ridge curvature. Their full-dimensional matrices are compared directly. Positive coefficients provide the AR1 eigenvalue lower bound and condition-number upper bound recorded in the object. These are labelled as bounds.

If `s_R=2psi_half-psi_full` and `g1=-X1's_R/m`, the centered score covariance is

\[
\Sigma_h=aa^{\mathsf T}\operatorname{Var}(g_1)
+C E(\|s_R\|_2^2/m^2).
\]

The half-bandwidth and full-bandwidth cross term is retained. The profile builder supports a random intercept and a random intercept with time slope, including a fixed positive-definite non-diagonal Lambda, under its stated conditional independence conditions.

## Profile target object

`build_profile_pair_conditional_scale_target_v2()` creates `asset_kind="paired_conditional_scale_target"` with `population_engine="conditional_gaussian_scale_v1"`. Each independent batch contains the declared number of genuine clusters. The existing scalar streamed profile solver computes roots at the target bandwidth, half that bandwidth, and the analysis bandwidth. The target bandwidth is `min(h_analysis/4,n_population^(-1/3))` unless explicitly supplied.

The scalar root is embedded as `beta0+delta e1` in every registered coordinate. Four full embedded estimates determine their mean, coordinatewise MCSE and between-batch variation. Coordinates beyond the first have their analytically identified target values. The object retains the target-to-refinement and target-to-analysis movements, conditional full score, nuisance diagnostics and curvature identities. Completed scalar roots and completed repeat records can be checkpointed independently.

`validate_profile_pair_conditional_scale_target_v2()` verifies the complete scenario identity and source dependencies. It reconstructs target embeddings, MCSE and bandwidth movements from actual repeat records. It also checks the full score identity and the two curvature coefficients against their ridge-identity counterparts. The numerical limits are the existing profile population limits. Final scale requires at least four independent batches of at least 100000 clusters each.

## The own marginal QR target

SQR uses the pooled observation-level QR population projection of the supplied centered fixed-effect columns. Its target satisfies

\[
E\left[X_1\left\{F_{b+\sigma(X_1)(U-q_\tau)\mid X_1}
(\delta_{QR}X_1)-\tau\right\}\right]=0.
\]

The Gaussian conditional mean identity embeds this root as `beta0+delta_QR e1`. `build_sqr_marginal_scale_target_v2()` evaluates the registered innovation CDF, including its probability clipping, against the Normal random-intercept distribution. Two nested quadrature precision levels check the same population root. This builder currently requires the Normal random-intercept law used in B16 and B19.

The resulting object has `asset_kind="sqr_marginal_qr_target"` and `population_engine="marginal_scale_quadrature_v1"`. The full target is available as `beta_target` and `beta_marginal`. Its stochastic target MCSE is zero because the target is computed by deterministic population integration; numerical root agreement and population equation residuals are recorded separately.

Four independent genuine-cluster batches audit the same marginal score. Their pooled-visit ratio uses independent-cluster MCSE, retaining variable visit counts and the shared random intercept. The score and interval target remain those of SQR. Experiment-level distances continue to use the declared profile target, while SQR coordinate interval coverage uses its own target.

## Source identity and storage

All new functions reside in simulation part08. Existing simulation parts01 through07 and the profile and native numerical sources retain their definitions. The original finite-X and structural conditional-Gaussian asset validators remain available for their object types.

Each new object records its exact numerical source dependencies and implementation commit. A descendant implementation may consume it when those dependencies are unchanged. Final construction requires committed source dependencies. The new paths are `construction_v4/targets/{id}.rds` for profile targets and `construction_v4/comparator_targets/SQR-DEBIASED-IID/{id}.rds` for SQR targets. Existing construction_v3 objects retain their paths and identities.
