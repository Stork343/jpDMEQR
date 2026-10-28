# Theory for Joint Profile Inference with Conditional Covariate Selection

## Parameter and population geometry

Let `ell_b(D,beta)` be the complete cluster profile loss in `METHOD.md`, `R_b(beta)=E ell_b(D,beta)`, and `beta_star=argmin R_0(beta)`. The target coordinate is `k`. The unsmoothed population Hessian is `H_0=R_0''(beta_star)`. Define

\[
 \omega_0=H_0^{-1}e_k,\quad
 c_0=\omega_0/\omega_{0,k},\quad
 J_0=c_0^TH_0c_0=1/\omega_{0,k}.
\]

For a selected set `A` containing `k`, write `beta_A^P` for the minimiser of `R_0` with coefficients outside `A` set to zero. This selected model parameter is used in the proof. The reported target remains `beta_star_k`.

The population density section theorem in the project source `research/profile-score-theory-20261022/main.tex`, labelled `thm:primitive-bandwidth`, gives

\[
 E g_b(\beta)=\nabla R_0(\beta)+bB_1(\beta)+O(b^2),
 \qquad E H_b(\beta)=R_0''(\beta)+bC_1(\beta)+O(b^2). \tag{1}
\]

The directional forms of these expansions follow by replacing design columns by linear contrasts. They are uniform in the sparse coefficient neighbourhood used below.

## Distribution and curvature conditions

Clusters are iid, with bounded sizes and fixed nuisance dimension one or two. `Lambda` has eigenvalues bounded above and away from zero. For dimension two the nuisance design contains an intercept and a bounded time coordinate, with bounded densities for pairwise time distances. The conditional outcome density has bounded density sections and their integrated derivatives through order three. Conditional unit outcome contrasts have bounded densities.

For every unit vector `v`, the conditional sub Gaussian norm of `X_j^T v`, given the nuisance design and cluster size, is bounded by a common constant. Sparse design eigenvalues are bounded above and away from zero. The complete population profile risk has positive restricted curvature and bounded upper curvature on a ball containing the considered sparse minimisers. Its first three directional derivatives satisfy the bounds supplied by the density sections and design moments. The risk is finite, and the average absolute response has a bounded expectation.

The population coefficient has sparsity at most `s`, with bounded Euclidean norm. The contrast `c_0` has sparsity at most `s_c` and bounded L1 norm. Its influence variance is positive and finite. The number of requested target coordinates is fixed. The two estimation folds and three inner selection folds have positive limiting sampling proportions.

Let `M` be the largest inner search budget, with the smallest budget of the same order, and put `D=7M+s+s_c+d_0+2`, `L=log(2np)`. The quantity `d_0` counts always included coordinates. The design regression adds at most `M` variables to the union of the three outcome and profile graph supports. The search budget satisfies `M/(s log n)->infinity` in the explicit sparse corollary, and admits a sufficiently large fixed multiple of `s_c` for graph approximation.

## Learning approximation

The search and graph proofs are given in the sections from sparse gradient concentration through learned contrasts in `THEORY.md`. Their objects before the final path construction are precisely the inner learners used here. The following argument states how those bounds apply to the actual selected sets, including the added design regression.

The cluster dual problem has bounded dimension. Its minimum norm optimal solution is piecewise affine in `Y-X beta`. On any support with at most `D` coordinates, its projected score subgraphs are described by a finite collection of bounded degree polynomial inequalities. Polynomial sign counting gives covering logarithm `O(D log(Cp/epsilon))`. The bounded dual box and conditional sub Gaussian design contrasts give a bounded directional second moment and an observed envelope of order `sqrt(DL)`. Symmetrisation, the entropy integral and Bernstein truncation yield

\[
 \sup_{\beta,v}|(P_T-P)v^Tg_0(\beta)|=O_p(\delta_n),
 \quad\delta_n=\sqrt{DL/n}+(DL)^{3/2}/n. \tag{2}
\]

The supremum is over the stated coefficient ball and unit sparse directions, including the selected spans. Nonunique dual faces can contribute through at most `D` independent clusters at a `D` dimensional coefficient under the joint outcome density. Their bound `CD sqrt(DL)/n` is included in (2). A normalised primal dual gap `epsilon_n` gives an epsilon subgradient. Evaluating its two subgradient inequalities at `beta +/- t v` and taking `t=sqrt(epsilon_n)` adds `C sqrt(epsilon_n)` to (2).

Set `nu_n=delta_n+n^(-1)+sqrt(epsilon_n)`. The selected span population gradient has norm `O_p(nu_n)`. Population curvature and the residualised design eigenvalue bounds imply that the largest eligible population gradient is at least `c sqrt(r_j/s)-C nu_n`, where `r_j` is the current excess risk. The computed dual ordering has that accuracy by (2). Bounded upper curvature then yields

\[
 r_{j+1}\leq(1-c/s)r_j+C\nu_n^2,
 \qquad r_M\leq e^{-cM/s}r_0+Cs\nu_n^2. \tag{3}
\]

For each trained path, the heldout loss difference is bounded by the average absolute design contrast. Its conditional sub Gaussian scale is proportional to the coefficient error. A finite union bound over the inner paths and population curvature give the validation oracle remainder `C log(CM)/n`. The chosen inner fits therefore have Euclidean error of order

\[
 e_n=e^{-cM/(2s)}+\sqrt{s}\nu_n+\sqrt{\log(CM)/n}. \tag{4}
\]

Their supports are retained in the union, so each final selected model contains a coefficient comparison vector attaining this risk. Restricted population minimisation cannot increase that risk. Thus `||beta_A^P-beta_star||_2=O_p(e_n)`. Adding the design regression support preserves the comparison vector. This argument also accounts for refitting a union rather than repeating a separately selected full data path.

The complete Schur identity is

\[
 v^TH_{i,b}v=\min_a\{\|W_{i,b}^{1/2}(X_iv-Z_ia)\|_2^2+a^T\Lambda a\}.
\]

It is positive. The cluster curvature moment proof labelled `thm:oracle` in the project source extends to arbitrary unit design contrasts and gives

\[
 E|v^TH_{i,b}(\beta)w|^r\leq r!C^r b^{1-r}\|v\|_2^r\|w\|_2^r.
\]

Bernstein's inequality on sparse coefficient nets and unit sphere nets gives the restricted operator rate

\[
 a_n=\sqrt{DL/(nh)}+DL/(nh). \tag{5}
\]

At the sparse population contrast, the nodewise empirical gradient has order `a_n+e_n+h`. A geometric penalty grid with bounded adjacent ratios contains a common candidate of this order on every inner fit. If the population penalty anchor is smaller than this order, the zero graph is already a valid comparison. Otherwise the anchors are comparable, and the lower ratio endpoint is below the required order. The Lasso basic inequality and sparse upper eigenvalues give prediction error `C s_c(a_n+e_n+h)^2` and support size `C s_c` for that candidate. Its support fits within the declared graph budget.

Positive quadratic refits preserve this comparison risk. A useful validation reference is the sparse contrast `c_0`. With `Delta=c-c_0` and `delta=e_n+h+n^(-1)`, restricted curvature and the population gradient defect give

\[
 E\{\ell_{graph}(c)-\ell_{graph}(c_0)\}
       \geq\kappa\|\Delta\|_2^2-C\delta\|\Delta\|_2.
\]

Adding a sufficiently large common multiple of `delta^2` makes this excess nonnegative and bounds `||Delta||_2^2` by a constant times that shifted excess. The cluster curvature moments then bound the variance of each heldout loss difference by `C/h` times the shifted excess. The reference loss and shift are common across the candidate ratios on each validation sample and cancel in their selection. Bernstein's inequality gives the validation oracle remainder `C log(CM log n)/(nh)`, together with the already included `delta^2` comparison error. Taking the union over all possible outcome path sizes permits their selected size to enter this graph calculation. The selected inner graph supports consequently contain contrasts approximating `c_0` with error at most

\[
 r_n=C\{(1+\sqrt{s_c})(a_n+e_n+h+n^{-1})
                  +\sqrt{\log(CM\log n)/(nh)}\}. \tag{6}
\]

Every final set `A` contains such a comparison contrast `v_A` with target entry one. The additional design support cannot remove it. Its size increases the sparse neighbourhood bound from `6M` to `7M`; it does not require a new population target or a new curvature condition.

## Selected model target bound

For each selected set, let `Delta_A=beta_A^P-beta_star`. Restricted stationarity gives `v_A^T R_0'(beta_A^P)=0`. Expand this identity at the full target:

\[
 0=J_0\Delta_{A,k}+(v_A-c_0)^TH_0\Delta_A+O(\|\Delta_A\|_2^2).
\]

The remainder bound follows from the population directional derivative condition and the bounded norm of `v_A`. By (4) and (6),

\[
 |\beta_{A,k}^P-\beta_k^\star|=O_p(e_n^2+r_ne_n). \tag{7}
\]

This is a target approximation statement for the selected model. It permits imperfect outcome selection and imposes no lower bound on a nonzero outcome coefficient.

## Exact preservation through conditional covariates

Suppose that, conditional on the cluster design and size, the covariate rows are independent centred Gaussian with covariance `Sigma`, and independent of the complete error and nuisance vectors. Let `Y=X beta_0+Zb+epsilon`. For any fixed selected set `A`, Gaussian conditional projection gives

\[
 X_{A^c}=X_A\Sigma_{AA}^{-1}\Sigma_{A,A^c}+E_A,
 \qquad E_A\ \hbox{independent of }X_A.
\]

Consequently the selected complete profile population coefficient, at every bandwidth, is

\[
 \beta_A^P=\beta_{0,A}+\Sigma_{AA}^{-1}\Sigma_{A,A^c}\beta_{0,A^c}. \tag{8}
\]

To verify (8), substitute the Gaussian projection into the outcome. At the displayed coefficient the residual vector, including its profiled nuisance minimiser, is independent of the centred selected design. Its expected complete profile gradient is zero. Positive population curvature gives uniqueness. The same calculation applies at bandwidth zero and at positive bandwidths.

If `supp(Sigma^(-1)e_k)` is contained in `A`, the target entry of the second term in (8) is zero. Indeed, solving the `A` block of `Sigma Sigma^(-1)e_k=e_k` gives `v_A=Sigma_AA^(-1)e_k`, while its complementary block gives `Sigma_Ac,A v_A=0`. Transposition proves `e_k^T Sigma_AA^(-1)Sigma_A,Ac=0`. Hence

\[
 \beta_{A,k}^P=\beta_{0,k}. \tag{9}
\]

This population result explains retaining conditionally associated covariates. The conditional neighbours can protect the target even when outcome prediction omits other nonzero coefficients. The simulation audit evaluates (8) from the declared Gaussian design after fitting; the generated coefficients do not enter selection or estimation.

## Paired profile estimation

Conditional on a fold's selection sample, its selected set is fixed and its estimation clusters retain their original iid law. Let `beta_A,b^P` be the selected population minimiser at bandwidth `b`. Applying (1) and the inverse function theorem on the selected coordinates gives

\[
 \beta_{A,b}^P=\beta_A^P-bH_{A,0}^{-1}B_{1,A}(\beta_A^P)+O(b^2).
\]

The expansion is uniform on the selected sets and coefficient ball. Using the same set at both bandwidths yields

\[
 2\beta_{A,h/2,k}^P-\beta_{A,h,k}^P
       =\beta_k^\star+O_p(h^2+e_n^2+r_ne_n). \tag{10}
\]

The prediction ridge in the implementation adds order `n^(-2)` to this expression. The empirical design eigenvalue condition bounds its gradient on the considered ball. Its contribution is negligible after multiplication by `sqrt(n)`.

On each estimation sample the two coefficient fits solve their entire selected score equations. Their influence rows are obtained from the corresponding selected Schur Hessians. Let `d` bound the selected dimension. Conditional concentration on coefficient and direction nets within that fixed selected space gives coefficient error of order `sqrt(d log n/n)` and curvature error of order `sqrt(d log n/(nh))+d log n/(nh)`. Their product is negligible on the scalar root scale when `d log n/sqrt(nh)->0`.

For a common sufficient condition covering fitted second moments, require `(d log n)^3/(nh)->0`. One way to verify it is to truncate the selected design norm at `C sqrt(d log n)`, apply the curvature moment bound to squared directional entries, and extend sparse coefficient and sphere nets. The variance term has size `sqrt(d log n/(nh^3))`; the bounded term has size `d^3(log n)^3/(nh^2)`. Under the displayed condition their sum is `o(1/h)`. Integrating the curvature along fitted coefficient segments then bounds score substitution in mean square by the squared coefficient error divided by `h`.

## Inference theorem

Under the distribution and learning conditions above, suppose

\[
 \sqrt n\{h^2+e_n^2+r_ne_n\}\to0,
 \quad e_n^2/h+r_n^2\to0,
 \quad (d\log n)^3/(nh)\to0. \tag{11}
\]

Use numerical tolerances whose score contribution is `o(n^(-1/2))`, as ensured by the decreasing tolerances in `small_profile.R`. Then

\[
 \sqrt n(\widehat\theta_k-\beta_k^\star)
  =-n^{-1/2}\sum_{i=1}^n\omega_0^Tg_{i,0}(\beta^\star)+o_p(1), \tag{12}
\]

and the variance computed by the algorithm satisfies

\[
 n\widehat s_k^2\ \longrightarrow_p\
 \sigma_k^2=\operatorname{Var}\{\omega_0^Tg_{i,0}(\beta^\star)\}.
\]

Therefore `(theta_hat-beta_star_k)/se -> N(0,1)` and the stated interval has limiting coverage `1-alpha`.

**Proof.** Conditional on the selection sample of one fold, apply the selected coefficient score expansion described above. Its target component equals the corresponding population coefficient plus the average selected influence contribution and an `o_p(n^(-1/2))` remainder. Population restricted curvature and the comparison contrast `v_A` imply that the selected inverse row approaches `omega_0` in Euclidean norm at order `r_n+e_n+h`.

To replace the fitted selected population scores by full target scores, centre each score at its own population expectation. Conditional mean square differences are bounded by `C(e_n^2/h+r_n^2+h)`, using the directional curvature moments, the density section limit and the bounded score moments. They tend to zero under (11). Conditional Chebyshev bounds control their root sample size scaled means. This conditioning is applied separately on each fold; it does not condition jointly on the two selected sets.

For each bandwidth, the leading centred influence contribution is thus the same independent cluster sum with row `omega_0`. Combining the two fits with weights `(2,-1)` retains that contribution once. Combining the folds with weights `n_f/n` gives each estimation cluster weight `1/n`. Equation (10) and the first part of (11) remove the population approximation remainder, giving (12). The independent cluster central limit theorem applies to the bounded dual score multiplied by a sub Gaussian design contrast.

For the variance, the recorded paired contribution differs in average squared norm from `-omega_0^T g_i,0(beta_star)` by `o_p(1)`. The fitted score substitution bound, convergence of the inverse rows and `e_n^2/h+r_n^2->0` establish this statement. Cauchy Schwarz controls the cross product in second moments. Centring within a fixed number of proportional folds preserves the limiting second moment. The iid law of large numbers proves the variance assertion, and Slutsky's theorem proves the confidence statement. The prediction ridge and certified numerical residuals enter only the already bounded vanishing remainder. \(\square\)

## Explicit high dimensional regime

For bounded `s`, `s_c` and `d_0`, the implemented asymptotic budget is `M=ceil(log(n)^3)`. With `h=n^(-3/10)` and the declared optimisation gaps, a sufficient regime is

\[
 \log(2np)=o\{n^{7/20}/(\log n)^3\}. \tag{13}
\]

In this regime `e_n=O(sqrt(DL/n))` and `r_n=O(sqrt(DL/(nh))+h)`. The largest target approximation product in (11) is controlled by `DL/sqrt(nh)->0`. The smoothing term is `sqrt(n)h^2=n^(-1/10)`. The estimation dimension satisfies `d=O(log(n)^3)`, so its second moment condition in (11) also holds. This gives a nonvacuous model class with exponentially growing dimension.

Growing sparsities retain the explicit bounds (4), (6) and conditions (11). Approximate sparse coefficients or contrasts enter through additional population comparison errors in those two bounds. The exact population statement (9) separately describes when the design regression protects a coefficient despite omitted outcome variables.
