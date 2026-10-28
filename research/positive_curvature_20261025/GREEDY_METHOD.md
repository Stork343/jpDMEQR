# Dual Guided Search and Cross Fitted Profile Score Inference

Method identifier: `greedy_crossfit_profile_score_20261025`.

The estimand is the unsmoothed regularised cluster profile coefficient with fixed `Lambda`. The algorithm searches the complete fixed effect design. Independent sampling units are clusters. The working nuisance design supports both an intercept and an intercept with a time slope.

## Search by complete design dual scores

On a training sample of `n_T` clusters, standardise each predictor by its cluster weighted RMS. Let `S` contain the always included coordinates. For each selected set solve

\[
 \min_{\theta_S,\gamma}
 \frac1{n_T}\sum_i\left\{\frac1{m_i}\sum_j
 \rho_\tau(Y_{ij}-X_{s,ij,S}^T\theta_S-Z_{ij}^T\gamma_i)
 +\tfrac12\gamma_i^T\Lambda\gamma_i\right\}
 +\frac{\|\theta_S\|_2^2}{2n_T}.
\]

Every omitted coefficient is zero for that search model. The deterministic initial ridge tends to zero with the number of training clusters. It makes each selected coordinate problem strongly convex and supplies a computable primal dual gap for that selected coordinate problem.

At the current dual vector `u`, compute `X_s^T u/n_T` across all columns. Add the unselected coordinate with largest absolute value, and jointly optimise the enlarged selected design and all nuisance effects. Warm starts retain the preceding coefficients and dual vector. The new coordinate starts at zero. The selected convex objectives cannot increase beyond their recorded optimisation error.

Five whole cluster folds construct these paths. The allowed numbers of added coordinates range from zero to
`min(p_penalised, floor(n_min_training/log(n_min_training))-p_always)`;
the implementation retains at least one addition when there are penalised variables. The number of variables is chosen by the smallest total exact heldout profile loss. The search never uses generated supports or coefficients. The full sample follows its own search path to that chosen size.

## A finite sample extension bound

Write `P_S` for the current selected primal objective multiplied by `n_T`, and let `D_S(u)` be its dual value. With zero L1 penalties,

\[
 D_S(u)=Y^Tu-\tfrac12\|X_{s,S}^Tu\|_2^2
 -\tfrac12\sum_i(Z_i^Tu_i)^T\Lambda^{-1}(Z_i^Tu_i).
\]

For any additional coordinate set `A` disjoint from `S`,

\[
 D_{S\cup A}(u)=D_S(u)-\tfrac12\sum_{j\in A}(X_{s,j}^Tu)^2.
\]

Every box feasible `u` remains dual feasible after adding coordinates because the ridge conjugate is finite. If the current normalised gap is `epsilon_S`, the improvement attainable by this extension is at most

\[
 \epsilon_S+\frac1{2n_T}\sum_{j\in A}(X_{s,j}^Tu)^2.
\]

For an extension of at most `r` coordinates, the sum of the `r` largest squared dual correlations gives an upper bound. The search orders additions by this certified upper bound on improvement. The argument applies directly to the nonsmooth profile objective.

## Heldout score and positive curvature

The selected heldout coefficient fits already exist in the search paths. For each cluster use the coefficient fitted without that cluster's fold. Evaluate its complete nuisance profiles at `h=n^{-3/10}` and `h/2` and combine their scores as `2*g_half-g_full`.

The positive single bandwidth Schur complement is represented by rectangular cluster factors. Nodewise L1 paths learn the requested coordinate directions. Their five cluster fold validation, objective normalisation and KKT checks follow L1_REFERENCE.md. The factor has N rows and p columns.

## Projected score inversion

For a target coordinate `k` and candidate value `t`, replace coordinate `k` in every heldout coefficient vector by `t`. Retain its other coefficients and the fitted direction `omega_k`. The projected cluster score is

\[
 U_{ik}(t)=\widehat\omega_k^T G_{i,h}^R
 (\widehat\beta_{-f,-k},\beta_k=t).
\]

Set

\[
 Z_k(t)=\frac{\sqrt n\,\overline U_k(t)}
 {\{n^{-1}\sum_i(U_{ik}(t)-\overline U_k(t))^2\}^{1/2}}.
\]

The point estimate is a local root of `Z_k(t)=0`. The reported interval is the connected component of `abs(Z_k(t))<=z_0.975` around that root. Outward searches and repeated interior refinements locate its two boundaries. The initial one step interval supplies only a starting width for numerical bracketing. Brackets expand deterministically; they do not use the generated target.

The algorithm checks a positive local projected derivative at the root and evaluates the acceptance component on interior grids. These are numerical diagnostics. The theorem requires local identification and a connected limiting acceptance region. Root residuals, derivatives, evaluation counts and the interior diagnostic grid are saved. The reported standard error is the local score standard error divided by the actual projected derivative at the root; confidence limits come from score inversion.

For one nuisance coordinate, exact profile evaluation sorts the scalar residual knots and solves the piecewise linear stationarity equation including ridge curvature. It reconstructs box feasible dual scores at active knots and checks the primal dual gap. For two nuisance coordinates, whitening the ridge reduces the quadratic term to a Euclidean norm. Restricting the criterion to each residual zero line gives an exact scalar problem. The algorithm examines those restricted optima and the stationary points in adjacent cells, with dual feasibility and gap checks for every accepted candidate. The finite search and its completeness argument are in FINITE_PROFILE.md.

Before root solving, compute `X*omega_k` and the predictor omitting coordinate `k`. Each subsequent score evaluation uses only low dimensional nuisance solves and scalar row projections. Its cost is independent of `p` once these products are available.

## Validation risk statement

Let `R(beta)` be excess population profile risk. Suppose heldout loss differences have conditional sub Gaussian scale bounded by `C*||beta-beta_star||_2`, and population curvature gives `R(beta)>=kappa*||beta-beta_star||_2^2/2` on all considered paths. For a fixed number `F` of folds and at most `M` candidates per path, a conditional concentration inequality on each fold and a union bound yield

\[
 |\widehat R_f(\beta_{f,m})-R(\beta_{f,m})|
 \leq \eta R(\beta_{f,m})+
 C_\eta\log(2FM/\delta)/n_f
\]

simultaneously for all folds and candidates, with probability at least `1-delta`. The empirical excess is computed relative to the same true coefficient on each validation fold; this reference cancels when selecting the candidate. Consequently

\[
 \sum_f\frac{n_f}{n}R(\beta_{f,\widehat m})
 \leq\frac{1+\eta}{1-\eta}
 \min_m\sum_f\frac{n_f}{n}R(\beta_{f,m})
 +C_\eta\frac{F\log(2FM/\delta)}n.
\]

The profile loss difference is bounded by the average absolute predictor difference because check loss is Lipschitz and the same nuisance penalty enters both criteria. A bounded cluster size and sub Gaussian fixed effect design therefore supply the stated envelope without requiring sub Gaussian responses. The curvature condition and a useful approximation along the search path remain model requirements.

If the selected models have at most `d_n` coordinates and the target has `s` nonzero coordinates, this risk bound implies an L1 error bound through `||Delta||_1<=sqrt(d_n+s)*||Delta||_2`. A primitive approximation rate for the dual guided path under the intended clustered model remains to be established.

## Conditional score inference

At `t=beta_star_k`, each fold's coefficient error has coordinate `k` equal to zero. The decomposition in THEORY.md applies to these nuisance only errors. Under its initial error, direction error, profile expansion and score substitution conditions, the centred leading term is the independent cluster sum with direction `H_0(beta_star)^{-1}e_k`. The remaining projected mean is bounded by
`mu_n*e_n+(h+c_n)*e_n+e_n^2+h^2`, plus the estimated direction stochastic contribution.

The displayed rate conditions in THEORY.md give `Z_k(beta_star_k) -> N(0,1)`. Test inversion gives pointwise coverage for the local acceptance interval under the identification conditions stated above. The primitive rate analysis for the actual search and direction selection rules is part of the remaining mathematical work. No coverage multiplier appears in the algorithm.

## Computational structure

For a path with `M` added coordinates, complete design scans cost order `M*N*p`. At step `d`, a primal dual iteration costs order `N*d+N*q+n*q^2`. The rectangular profile factor and requested nodewise paths require order `N*p` storage. The number of path iterations and scalar root evaluations is measured in every run.

## Related methods

Greedy selection, primal dual linesearch, nodewise inversion and cross fitting have established precedents. Ing and Lai study OGA in sparse linear regression. Cha, Chiang and Sasaki combine greedy estimation with debiased inference. Malitsky and Pock supply the optimisation algorithm used for each convex subproblem. The proposed research concerns their development for the complete bounded cluster quantile profile, including exact dual bounds, profile smoothing, positive curvature and cluster score inversion. A literature comparison is required before claiming priority.

- [Ing and Lai 2011](https://www3.stat.sinica.edu.tw/statistica/oldpdf/a21n41.pdf)
- [Cha Chiang and Sasaki](https://arxiv.org/abs/2108.09520)
- [Malitsky and Pock 2018](https://arxiv.org/abs/1608.08883)
- [Van de Geer and colleagues 2014](https://arxiv.org/abs/1303.0518)
- [Chernozhukov and colleagues 2018](https://academic.oup.com/ectj/article/21/1/C1/5056401)
