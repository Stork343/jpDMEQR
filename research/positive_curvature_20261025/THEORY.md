# Optimisation and Conditional Inference

All statistical formulas in this document use the original coefficient units. RMS scaling in the implementation is an invertible change of variables. With coefficient scale matrix `D`, the raw direction equals `D^{-1} omega_scaled / D_kk`; this transformation exactly preserves the coordinate correction and its variance. Random empirical scales are required to be bounded above and away from zero on the model class.

## Joint convex optimisation

Let `P_T` be the training objective in L1_REFERENCE.md multiplied by its number of clusters `n_T`. Setting its L1 penalties to zero gives the selected orthogonal coordinate problem in ORTHOGONAL_METHOD.md. The identity `rho_tau(r)=max_{s in [tau-1,tau]} s*r` gives

\[
 P_T(\theta,\gamma)=G(\theta,\gamma)
 +\max_{u\in C}\{Y^Tu-u^TA(\theta,\gamma)\},
\]

where `A(theta,gamma)=X_s theta+Z_block gamma` and

\[
 G=\tfrac12\|\theta\|_2^2+\sum_j n_Ta_j|\theta_j|
            +\tfrac12\sum_i\gamma_i^T\Lambda\gamma_i.
\]

The conjugate of the fixed effect summand is `S(v_j,n_T*a_j)^2/2`. The conjugate of the nuisance quadratic is `z_i^T Lambda^{-1} z_i/2`. Fenchel duality therefore gives the displayed dual objective in L1_REFERENCE.md. The primal is coercive because its quadratic terms are positive definite. The dual box is compact and contains zero. A saddle point exists, and the primal minimiser is unique.

Algorithm 1 of Malitsky and Pock applies with `K=-A` and `f*(u)=indicator_C(u)-Y^T u`. Both proximal maps are evaluated exactly. Its linesearch condition is checked on the concatenated fixed effect and nuisance adjoints. Steps are bounded above by ten. The finite valued continuous function `G` meets the continuity condition of their convergence theorem, so the iterates converge to the unique primal solution and the computed gap tends to zero. The iteration cap is a numerical termination rule; reaching it raises an error.

If the normalised primal dual gap is `epsilon_T`, the returned coefficient vector is an `epsilon_T` approximate minimiser of the profiled penalised objective. To see this, minimising over nuisance coefficients can only decrease the returned primal objective, and the dual value is a lower bound on the joint optimum. Thus no inner nuisance tolerance needs to be accumulated across initial iterations.

Under a weighted restricted Bregman curvature condition, bounded positive penalty loadings, and score domination including the initial `1/n_T` ridge gradient, the approximate Lasso basic inequality gives

\[
 \|\widehat\beta-\beta^\star\|_1
 \lesssim s\lambda+\sqrt{s\epsilon_T}+\epsilon_T/\lambda.
\]

The inequality uses an augmented cone with slack proportional to `epsilon_T/lambda`. Curvature must hold on that cone and on a radius containing the optimiser. Always included coordinates enter the augmented support. This L1 bound applies to the reference initial estimator under its score domination and curvature premises. The orthogonal search uses the validation risk inequality in GREEDY_METHOD.md. Its primitive path approximation rate remains a research requirement.

## Positive curvature factor

For any vector `v`, the Schur identity yields

\[
 v^TX_i^TL_iX_iv
 =\min_a\{\|W_i^{1/2}(X_iv-Z_ia)\|_2^2+a^T\Lambda a\}\ge0.
\]

Thus `L_i` is positive semidefinite, and its symmetric square root exists, including when some kernel weights are zero. `F_i=L_i^{1/2}X_i` satisfies `F_i^TF_i=H_{i,h}`. The implementation removes eigenvalues of the small `m_i` matrix only within its declared roundoff tolerance and raises an error for a larger negative value. No change to the statistical curvature is made.

## Nodewise identities

Work in the unit diagonal curvature coordinates, with sample Gram matrix `C`. For target `k`, let `a` minimise `||F_k-F_-k a||^2/(2n)+lambda*||a||_1`. Its KKT equations give `C_-k,-k a-C_-k,k+lambda*z=0` with `z` a valid L1 subgradient. Put `c_k=1`, `c_-k=-a` and

\[
 t^2=c^TCc+\lambda\|a\|_1=C_{kk}-C_{k,-k}a.
\]

Consequently `(C c/t^2)_k=1` and the off diagonal inverse defects are bounded by `lambda/t^2`. Rescaling gives the formula in L1_REFERENCE.md. A numerical KKT error adds its measured contribution to these identities. For a nonzero response column and positive penalty, the residual plus penalty expression for `t^2` is positive. Each selected direction and its actual inverse defect are recorded.

## Cross fitted decomposition

Let `I_1,...,I_F` be the fixed cluster folds, with `F=5` and sizes proportional to `n`. The number of initial ratio candidates is fixed at eight. For a selected candidate write `beta_f` for the estimate trained without `I_f`, `Delta_f=beta_f-beta_star`, and `Delta_bar=sum_f(n_f/n)*Delta_f`.

Define the paired score `G_i^R(beta)=2g_i,h/2(beta)-g_i,h(beta)`. The positive curvature estimate uses `H_i,h(beta_f)` on `i in I_f`. Let `omega_hat` be the nodewise direction and let `H_hat` denote that positive empirical curvature. Direct substitution gives the exact finite sample identity

\[
 \widetilde\beta_k-\beta_k^\star
 =-\widehat\omega^TG_n^R(\beta^\star)
 +(e_k-\widehat H\widehat\omega)^T\overline\Delta
 +\widehat\omega^T\left\{\widehat H\overline\Delta
 -\frac1n\sum_f\sum_{i\in I_f}
 [G_i^R(\beta_f)-G_i^R(\beta^\star)]\right\}.
\]

Assume the population profile expansion and local Lipschitz curvature established for the original unsmoothed profile target. Uniformly on the initial neighbourhood,
`H_h(beta)=H_0(beta)+O(h)` and `E G_i^R(beta_star)=O(h^2)`.
Assume the empirical single and paired curvatures differ from their population versions by `O_p(c_n)` on every fold and candidate neighbourhood. Suppose `max_f ||Delta_f||_1=O_p(e_n)` and `||omega_hat||_1=O_p(1)`.

Expanding each fold score along its own segment and comparing every population curvature to `H_0(beta_star)` bounds the final brace, after projection, by

\[
 O_p\{(h+c_n)e_n+e_n^2\}.
\]

This argument also accounts for using an average of fold coefficients with an average curvature: insert the common matrix `H_0(beta_star)` in both quantities. Its products cancel exactly; the remaining differences have the stated order. The positive curvature need not equal the derivative of the paired score in a finite sample.

## Conditional normality theorem

Let `omega_0=H_0(beta_star)^{-1}e_k`, with bounded L1 norm. Assume

\[
 \|\widehat H\widehat\omega-e_k\|_\infty=O_p(\mu_n),
 \quad \|\widehat\omega-\omega_0\|_1=O_p(r_n),
 \quad \|\widehat\omega\|_1=O_p(1).
\]

Require all preceding initial and empirical conditions uniformly over the eligible finite ratio family. Suppose the centred cluster score maximum is `O_p(sqrt(log(2p)/n))` and

\[
 \sqrt n\{\mu_n e_n+(h+c_n)e_n+e_n^2+h^2\}\to0,
 \qquad r_n\sqrt{\log(2p)}\to0.
\]

Assume the independent cluster influence contributions
`xi_i,h=-omega_0^T{G_i^R(beta_star)-E G_i^R(beta_star)}`
satisfy a triangular array Lindeberg condition and their variance tends to a positive finite `sigma_k^2`. Then

\[
 \sqrt n(\widetilde\beta_k-\beta_k^\star)
 =n^{-1/2}\sum_i\xi_{i,h}+o_p(1)
 \ \,\Longrightarrow\,\ N(0,\sigma_k^2).
\]

Proof. The exact decomposition bounds the inverse defect term by `mu_n*e_n` and the curvature remainder by `(h+c_n)*e_n+e_n^2`. Decompose the leading empirical score into its mean and centred part. Its mean contributes `O_p(h^2)`. Replacing the estimated direction by `omega_0` in the centred part contributes `O_p(r_n*sqrt(log(2p)/n))`. The displayed conditions make every scaled remainder vanish, leaving the independent cluster central limit theorem. Uniform bounds over the fixed finite family permit choosing the ratio by the observed validation losses.

## Variance consistency

Assume the uniform score substitution bound gives

\[
 \frac1n\sum_f\sum_{i\in I_f}
 [v^T\{G_i^R(\beta_f)-G_i^R(\beta^\star)\}]^2=O_p(e_n^2/h)
\]

for the L1 bounded directions used in the argument. Require `e_n^2/h -> 0`, `r_n -> 0`, a bounded maximum empirical coordinate second moment, and a law of large numbers for the true projected squared scores. The fitted projected scores then differ from the true ones in empirical mean square by `o_p(1)`. Cauchy Schwarz controls the cross term in their empirical second moments. Centring preserves convergence, giving consistency of the cluster sandwich and studentised normality.

## What is established by these arguments

The duality gap, positive factor identity, nodewise KKT identity and cross fitted decomposition are finite sample algebraic statements. The optimisation convergence follows from the cited primal dual theorem. Asymptotic normality and sandwich consistency follow under the displayed initial and direction rate conditions. Establishing those rates for the implemented model size and direction selection rules, over the intended model class, remains a research requirement.

## Projected score confidence sets

For candidate `t`, replace coordinate `k` in each heldout coefficient vector by `t`, retaining the other coordinates and the estimated direction. At the true value `t_0=beta_star_k`, the fold coefficient errors have coordinate `k` equal to zero. Insert `H_0(beta_star)` into the score expansion. Since `omega_0^T H_0(beta_star)=e_k^T`, the population linear term in these nuisance errors vanishes.

Under the preceding uniform conditions, the remaining mean and empirical terms have the orders in the conditional theorem. The projected, studentised score therefore satisfies `Z_k(t_0) -> N(0,1)`. Inverting `abs(Z_k(t))<=z_(1-alpha/2)` gives pointwise asymptotic coverage. If the population projected derivative is bounded away from zero in the identified neighbourhood and the empirical derivative and variance converge uniformly there, the consistent local root has the same leading cluster influence representation. Its local acceptance component has the stated limiting coverage.

The estimated direction is fixed while inverting a coordinate score. Multiplying that direction by any positive constant leaves the studentised score, root and acceptance set unchanged. This also leaves the local derivative based standard error unchanged because its numerator and denominator scale together.

The implementation computes the nodewise normalisation as `F_k^T r_k/n`. Its difference from `||r_k||^2/n+lambda*||a_k||_1` equals the fitted coefficient's inner product with its KKT residual. This identity is checked numerically and retains the optimisation residual in the direction diagnostics.
