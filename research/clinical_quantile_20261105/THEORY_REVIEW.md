# Review of the conditional quantile construction

## Target identification

The parameter is the coefficient of the stated conditional quantile model. Averaging the check loss within patients preserves this target under that model. Under approximation, the point estimator targets the explicitly defined patient weighted projection. The correlation assisted score has a conditional moment requirement beyond the projection definition. Both statements are necessary because covariance weighting can change a projection target under misspecification.

## Orthogonality calculation

Let a_i=C_i^(-1)F_i r_i and let Delta_i=d_i delta_theta+X_i delta_beta. Conditional expectation of the score at a perturbed quantile satisfies

\[
E[\widehat a_i^T u_i(q_i+\Delta_i)\mid W_i,m_i]
=-\widehat a_i^TF_i\Delta_i+R_i,
\quad |R_i|\leq L\sum_j |\widehat a_{ij}|\Delta_{ij}^2/2.
\]

The population projection makes E[a_i^T F_i X_i]=0. Therefore the linear nuisance term is E[(estimated a_i-a_i)^T F_i X_i delta_beta]. Conditional Cauchy Schwarz bounds it by the product of instrument and nuisance prediction errors when the density and instrument moments are controlled. The theta term is cancelled by the one step update. Its remaining contribution is the preliminary theta error multiplied by the sensitivity estimation error.

Cross fitting makes the score evaluation patients independent of the learned functions for their fold. If the estimated influence function converges in L2, the conditional variance of its centred empirical difference is o(1/n). A fixed number of folds preserves this order by the triangle inequality, even though their training samples overlap. The population remainders must separately be o(n^(-1/2)); cross fitting does not eliminate them.

This verifies the expansion conditional on the rates in METHOD.md. A primitive rate theorem for the adaptive elastic net and density spacing fits remains to be established. In particular, prediction consistency alone is insufficient for the normal approximation.

## Covariance calculation

At the true conditional density and score covariance, E[phi_i^2]=E[r_i^T F_i C_i^(-1)F_i r_i]=J. The resulting influence variance is 1/J. For a working covariance, the variance is E[phi_i^2]/J^2 and is estimated from complete patient contributions. Observations within a patient must remain combined before computing variance.

The numerical checks compare the analytic exchangeable transformation with direct matrix inversion on actual GNHS designs. They also verify the projection derivative and the equality between sensitivity and variance in the density constant oracle calculation. These are algebraic checks, not tests of clinical coverage.

## Computational certificate

For lambda>0 and 0<alpha<1, the conjugate of lambda alpha |b|+lambda(1-alpha)b^2/2 is (|g|-lambda alpha)_+^2/{2 lambda(1-alpha)}. The restricted fit dual must obey the full coordinate conditions, including omitted variables. Adding every violating coordinate in finite batches terminates after at most all candidate coordinates have entered. The certificate is evaluated on the full design, so a sparse active computation does not impose a statistical support cap.

The alpha=0 and alpha=1 endpoints have their own quadratic and linear optimisation formulations. They optimise the same patient weighted loss with the corresponding declared penalties. Exact selected fits are checked with Clarabel on the full design.

## Required scientific comparisons

The point estimate comparison must use the same patients and inputs as quantile Lasso, mean elastic net calibration and clinical quantile regression. The inference comparison holds the preliminary quantile fit fixed and changes only the working covariance and its associated projection. Nominal coverage and interval length need dedicated experiments with known conditional quantile coefficients. The empirical prediction runs cannot establish that coverage result.
