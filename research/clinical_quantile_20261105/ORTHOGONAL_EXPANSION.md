# Patient cross fitted quantile inference

## Model and normalisation

At sample size n, observe n independent identically distributed patients. A patient has at most M observations, where M is fixed. Conditional on its complete design and measurement count, each outcome satisfies

\[
P(Y_{ij}\le d_{ij}\theta_0+X_{ij}^T\beta_0\mid W_i,m_i)=\tau.
\]

The conditional densities are bounded, have a uniformly bounded derivative on a common neighbourhood of the target quantiles, and are positive there. The quantities below may depend on n through the dimension and the distribution. Let F_i contain these densities. For a positive definite working covariance, define the projection and sensitivity in METHOD.md and write the normalised instrument as a_i=C_i^(-1)F_i r_i/J. It satisfies

\[
E[a_i^TF_i d_i]=1,\qquad E[a_i^TF_iX_i]=0.
\]

Let K be fixed and let every evaluation fold contain a proportion bounded away from zero of the patients. Its complementary training patients determine the complete feature transformations and nuisance fits. The estimates and instruments are represented in the original outcome and target variable units before aggregation.

## Conditions on learned quantities

Conditional on the training patients for each fold, suppose the following bounds hold with probability tending to one. Norms refer to a fresh patient from the same population.

* The preliminary target error is O(a_n), and the L4 norm of X_i times the nuisance coefficient error is O(a_n). The target design has a bounded fourth moment, so the complete quantile prediction error also has L4 norm O(a_n).
* The L2 error of the normalised instrument is O(b_n), and the learned instrument has a bounded fourth moment.
* Both a_n and b_n tend to zero, and sqrt(n){a_n squared+a_n b_n} tends to zero.
* The oracle influence variance is bounded above and away from zero. Its fourth moment is bounded uniformly in n.

These conditions include the error of the estimated sensitivity through the normalised instrument. A consistent, positive sensitivity together with an unnormalised instrument error rate and a sensitivity error rate gives the corresponding b_n by the triangle inequality.

Define u_i(q) componentwise as tau-I(Y_ij<=q_ij). For fold f let q_hat_f=d_i theta_hat_f+X_i beta_hat_f and let a_hat_f be its normalised instrument. The estimator is

\[
\widehat\theta=\sum_f\frac{n_f}{n}
\left\{\widehat\theta_f+\frac1{n_f}\sum_{i\in I_f}\widehat a_{if}^T u_i(\widehat q_f)\right\}.
\]

## Expansion

Under the stated conditions,

\[
\sqrt n(\widehat\theta-\theta_0)
=n^{-1/2}\sum_{i=1}^n a_i^Tu_i(q_0)+o_p(1).
\]

To verify the population remainder, set Delta_i=d_i delta_theta+X_i delta_beta for a fold. Conditional expectation of the estimated score has the expansion

\[
E[\widehat a_i^T u_i(q_0+\Delta_i)]
=-E[\widehat a_i^TF_i d_i]\delta_\theta
-E[\widehat a_i^TF_iX_i\delta_\beta]+R_f.
\]

The density derivative bound and conditional Cauchy Schwarz give |R_f|<=C E[||a_hat_i|| ||Delta_i|| squared]=O(a_n squared). After adding delta_theta from the preliminary estimate, the target term is

\[
-\delta_\theta E[(\widehat a_i-a_i)^TF_i d_i],
\]

which is O(a_n b_n). Orthogonality gives

\[
E[\widehat a_i^TF_iX_i\delta_\beta]
=E[(\widehat a_i-a_i)^TF_iX_i\delta_\beta]
=O(a_n b_n).
\]

Thus the conditional population remainder is o(n^(-1/2)). This step controls all three sources of first order sensitivity: the target starting value, nuisance regression and learned instrument normalisation.

For the empirical term, decompose the difference between the estimated and oracle influence into an instrument error times u_i(q_0), and a_hat_i times the change in the indicator vector. The first squared L2 norm is O(b_n squared). The density bound implies that the conditional probability of a changed indicator is at most C |Delta_ij|. The bounded fourth moment of a_hat_i and Cauchy Schwarz then bound the second squared L2 norm by O(a_n). Consequently the centred empirical difference on an evaluation fold has conditional variance o(1/n_f).

Conditional Chebyshev bounds this difference by o_p(n^(-1/2)). Combining a fixed number of folds preserves that order. Independence of the different training samples is unnecessary for this last step: each evaluation fold is independent of its own training sample, and a finite sum of terms that are o_p(n^(-1/2)) has the same order.

The patient central limit theorem applies to the oracle sum. The same L2 approximation and the fourth moment condition give consistency of the sample influence second moment. Subtracting the fold mean removes the fold specific conditional mean, which converges to zero. Therefore the standard error defined in METHOD.md consistently estimates the scale in this expansion.

## Two working covariances

Apply the preceding argument jointly to independence and exchangeable working covariances. The Cramer Wold argument gives a centred joint expansion with influence covariance S, and the paired patient covariance estimate is consistent.

If D=S_00+S_11-2S_01 is positive, the clipped scalar minimiser t_hat converges to t_star. Writing the combined estimate around t_star gives

\[
\widehat\theta_A-\theta_0
=(1-t_\star)(\widehat\theta_0-\theta_0)
+t_\star(\widehat\theta_1-\theta_0)
+(\widehat t-t_\star)(\widehat\theta_1-\widehat\theta_0).
\]

The final product is o_p(n^(-1/2)). The limiting variance is the minimum of (1-t,t)S(1-t,t)^T over t in [0,1], so it is no greater than either endpoint variance. When D=0, the component influence functions agree in L2; their estimate difference is o_p(n^(-1/2)), and every bounded convex combination has the same first order expansion.

## Preliminary smoothing and remaining rates

The centred kernel in CENTERED_SMOOTHING.md has zero first moment. Under positive local curvature and the stated density smoothness, its population prediction contribution is O(delta squared). With delta=n^(-1/5), that contribution is O(n^(-2/5)), and its squared remainder satisfies the displayed expansion condition. The registered numerical ridge contributes an order n^(-1) objective term for bounded comparison coefficients in its declared coordinates.

A full result for the implemented procedure still requires rates for the data selected regression, the quantile spacing density estimate and the weighted projection. The theorem above identifies exactly which norms and products those results must control. It does not deduce those rates solely from successful numerical convergence.
