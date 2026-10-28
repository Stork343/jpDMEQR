# Profile bandwidth pair

The target is the minimiser of the unsmoothed complete ridge-profile population risk, with the working matrix `Lambda` fixed. The mathematical development is in [`main.tex`](../../research/profile-score-theory-20261022/main.tex).

## Profile components

At a common coefficient vector and on the same clusters, solve the complete nuisance profiles separately at bandwidths `h/2` and `h`. Define

\[
G_{i,h}^{R}(\beta)=2g_{i,h/2}(\beta)-g_{i,h}(\beta),
\qquad
J_{i,h}^{R}(\beta)=2H_{i,h/2}(\beta)-H_{i,h}(\beta).
\]

Each component uses the original fixed-effect design in its gradient and the exact Schur-complement Hessian. These quantities differentiate the sample criterion `2 Q_(h/2) - Q_h`.

If the population gradient has the uniform expansion

\[
E g_{i,h}(\beta)=\nabla R_0(\beta)+h b_1(\beta)+O(h^2),
\]

then `E G_R(beta_star) = O(h^2)`. The mathematical development proves this expansion by partitioning the joint outcome integral according to the interior coordinates of the dual score. It assumes uniformly integrable density sections through the second derivative, bounded nuisance designs and a compact kernel positive inside its support. Independent Normal or standardised `t3` errors with location cluster effects and bounded positive scales satisfy the stated density conditions. The result applies to bounded clusters and fixed nuisance dimension.

For random intercepts and random time slopes, the same development establishes uniform Hessian concentration on sparse parameter neighbourhoods. Its squared-Hessian bound gives mean-square score plug-in of order `e_n^2/h` for a starting estimator with L1 error `e_n`. The coordinate theorem retains the initial estimator's support and error conditions, local information conditions, and feasible sparse precision rows. These quantities specify how a selected refit enters inference.

## Computational interface

`profile_bandwidth_pair_v2()` takes the data, a common `beta`, `tau`, `h`, and `lambda_gamma`. It returns the combined objective, gradient, Hessian and cluster scores, together with both constituent fits. It also combines the cluster-fold Hessian sums using the same weights. The score definition is recorded as `profile_bandwidth_pair_2_half_minus_full`.

Both nuisance profiles must converge at their declared gradient tolerance. A failed profile terminates evaluation and reports the two maximum nuisance gradients. Named fixed-effect columns identify the output coordinates. `hessian_min_eigenvalue` records the realised curvature. The combined Hessian is symmetric and can be indefinite; local curvature assumptions apply directly to this matrix.

`profile_bandwidth_pair_wald_v2()` accepts a paired component object, the common starting coefficient vector, explicitly supplied precision columns and their coordinate indices. A supplied precision column corresponds to the combined Hessian. The function reports its inverse defect and computes

\[
\widetilde\beta_k=\bar\beta_k-
\widehat\omega_k^{\mathsf T}\overline G_h^R,
\qquad
\widehat{\operatorname{se}}_k^2=
\frac1{n^2}\sum_i
\left[\widehat\omega_k^{\mathsf T}
(G_{i,h}^R-\overline G_h^R)\right]^2.
\]

The projected score is combined before centring. Its variance contains

\[
4\operatorname{Var}(g_{h/2})+\operatorname{Var}(g_h)
-4\operatorname{Cov}(g_{h/2},g_h).
\]

The function requires positive projected-score variances and finite inputs. Its interval interpretation uses the initial-estimator, precision, score-bias and plug-in conditions stated in the mathematical development.

## Reproduction

The [research entry points](../scripts/research/README.md) compute the analytic location example by numerical integration, evaluate the existing cluster data-generating function, check the matched differential and covariance relations, and generate the numerical tables in the mathematical development. Joint-density chart integration checks total probability and population curvature for two independent Normal observations. The tests cover both random intercepts and random intercepts with slopes, including a non-diagonal working penalty matrix.

## Fitting and inference

`fit_profile_pair_v2()` obtains a convex weighted L1 estimate at `h_est`,
including both the stochastic score contribution and the order-`h_est`
population-score contribution in its penalty. It refits the selected support
together with all unpenalised and prespecified target coordinates at `h_inf`.
The dimension limit is `floor(n/log(max(n,3)))`. The returned paired scores
and Jacobian are evaluated at this common refit coefficient vector.

`infer_profile_pair_v2()` selects each Dantzig tolerance by held-out inverse
defect on deterministic cluster folds. Its admissible domain is `0 < mu < 1`.
The requested grid is retained, with values at least one recorded as
`excluded_zero_direction`, because they admit the zero precision direction.
Inference uses the selected full-design direction and the paired cluster
score for both the correction and its uncertainty.

Registry rows declare `method_version=profile_pair_20261022` to use this
procedure. The practical route and the true-support, true-nuisance,
population-direction and reduced-inverse diagnostics retain the same paired
score definition. Population target and direction assets have separate
versioned paths and validators. Development assets and the small registry
checks do not qualify as final Monte Carlo outputs.

For a prespecified low-dimensional application design,
`lowdim_profile_pair_inference_v2()` uses an unpenalised convex refit followed
by the paired correction with the full reduced Jacobian inverse. It requires
positive identified reduced curvature. Its subject jackknife repeats this
entire estimator for every omitted subject at the fixed full-sample bandwidth.
The sandwich follows the stated paired-score theorem. Jackknife intervals are
reported as a prespecified sensitivity analysis; their variance consistency
also needs a uniform leave-one expansion.
