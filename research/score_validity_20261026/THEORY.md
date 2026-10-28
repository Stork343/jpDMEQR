# Statistical Theory for Honest Paired Roots

## Population parameter and notation

For an independent cluster `D=(Y,X,Z,m)`, define

\[
 \ell_b(D,\beta)=\min_\gamma\left\{
 m^{-1}\sum_j\rho_{\tau,b}(Y_j-X_j^T\beta-Z_j^T\gamma)
 +\tfrac12\gamma^T\Lambda\gamma\right\},\quad
 R_b(\beta)=E\ell_b(D,\beta).
\]

The target `beta_star` minimises `R_0`. Write `g_b` for the complete profile score and `H_b` for its Schur curvature. Set `H_0=R_0''(beta_star)` and let `c_0` be the population contrast with target entry one satisfying `(H_0 c_0)_{-k}=0`. Its directional curvature is `J_0=c_0^T H_0 c_0`. All rates use the number of clusters.

The definitions of `g_b` and `H_b` include the nuisance ridge. In particular, `g_b=-X^T psi_b/m`, and the Schur identity includes `A^T Lambda A`. The density section expansion and the cluster curvature moment argument are in the project source `research/profile-score-theory-20261022/main.tex`, Theorems labelled `thm:primitive-bandwidth` and `thm:oracle`. The arguments below use their directional forms, obtained by replacing design columns by linear combinations of columns.

## Structural conditions

Clusters are iid, `m<=m_max`, and the eigenvalues of the fixed matrix `Lambda` are bounded above and away from zero. The nuisance design is an intercept or an intercept with a bounded time coordinate. In the latter case, pairwise time distances have bounded densities. Conditional on the design, outcome vectors have a joint density whose sections and their derivatives through order three have uniformly bounded integrals, as defined in the cited project theorem. Conditional densities of unit outcome contrasts are uniformly bounded. Normal and standardised Student t errors with three degrees of freedom, mixed over bounded scale changes, satisfy the stated density requirements.

Conditional on `(Z,m)`, all linear design contrasts `X_j^T v` have sub Gaussian norm at most `C||v||_2`. The population design covariance has bounded sparse eigenvalues. The complete profile risk has restricted curvature between positive constants on a fixed Euclidean ball containing its considered sparse minimisers and search paths. This is a condition on population profile risk, including the nuisance penalty. Its derivatives in any three unit design directions are bounded. The density section argument and the design contrast moment bound give this derivative bound and the directional bandwidth expansions. The population risk is finite and `E m^{-1}sum_j|Y_j|` is bounded.

The target has at most `s` nonzero coordinates, bounded Euclidean norm, and lies inside the curvature ball. The contrast `c_0` has at most `s_c` nonzero coordinates and bounded L1 norm. Its influence variance `Var(c_0^T g_0)/J_0^2` is bounded above and away from zero. Population curvature at the target is invertible on the coordinates defining this contrast. Local curvature for the population nodewise regressions is bounded below. The conditions permit correlated designs and weak coefficients; no minimum nonzero coefficient is required.

Every fold has a fixed positive sampling proportion. The number of outer and inner folds is fixed. Let `M` be the smallest inner search budget, `D=6M+s+s_c+d_0+2`, where `d_0` counts always included coordinates, and `L=log(2np)`. The budget is large enough to contain the sparse comparison supports and obeys `M/(s log n) -> infinity` for the explicit rate corollary. A fixed large multiple of `s_c` fits within the graph budget.

## Sparse gradient concentration

On a fixed support with at most `D` coordinates, the unsmoothed cluster dual problem is a box constrained quadratic problem of dimension at most `m_max`. Its minimum norm optimal dual is piecewise affine in the residual vector. A description can be constructed by enumerating lower, upper and free dual coordinates, solving the finite KKT systems, and selecting their minimum norm solution. The number of pieces depends on `m_max`, rather than `n` or `p`. Singular faces can be subdivided using their ranks, with the same finite bound.

Consequently the class

\[
 \{v^Tg_0(D,\beta):\|v\|_2=1,
       |\operatorname{supp}(v)\cup\operatorname{supp}(\beta)|\leq D\}
\]

has polynomial subgraph growth on each support: evaluating a subgraph on `N` observations requires signs of `O(N)` polynomials of bounded degree in `O(D)` real parameters. The polynomial sign count is at most `(CN)^(CD)`. Taking the union over supports adds `CD log(2p)` to its logarithm. This gives covering logarithm `CD log(Cp/epsilon)` after normalising by its envelope. The polynomial sign bound includes equality cases.

The score envelope in a unit direction is bounded by `C max_j|X_j^T v|`. Its second moment is bounded, while the envelope over all eligible directions and observed clusters is `O_p(sqrt(DL))`. Symmetrisation, the entropy integral and Bernstein truncation therefore give

\[
 \sup_{\beta,v}|(P_n-P)v^Tg_0(\beta)|
 =O_p(\delta_n),\qquad
 \delta_n=\sqrt{DL/n}+(DL)^{3/2}/n. \tag{1}
\]

For completeness, at entropy level `u`, the variance term is `C sqrt(u/n)` and the truncated envelope term is `C sqrt(DL)u/n`; integrating the polynomial covering bound and taking `u` of order `DL` gives (1). Sub Gaussian tails make the truncation remainder smaller than these terms. Under a conditional joint outcome density, more than `D` independent clusters cannot simultaneously lie on nonunique dual faces at a `D` dimensional coefficient: each such face imposes an independent nontrivial affine outcome restriction after eliminating the cluster nuisance. The exceptional contribution for other optimal dual selections is at most `CD sqrt(DL)/n`, which is included in (1).

An approximate dual solution also has a direct bound. A normalised primal dual gap at most `epsilon` supplies an `epsilon` subgradient of the empirical profile risk at its returned coefficient. Evaluate its subgradient inequality at `beta+t v` and `beta-t v`. Integrating (1) along these two segments and using the bounded population Hessian bounds the difference from the population directional gradient by

\[
 C\{\delta_n+t+\epsilon/t\}.
\]

Taking `t=sqrt(epsilon)` yields `C(delta_n+sqrt(epsilon))`. This calculation accounts for the dual vector used to choose the next coordinate.

## Greedy approximation and inner risk selection

Let `r_j=R_0(beta_j)-R_0(beta_star)` at the `j`th restricted fit. Put

\[
 \nu_n=\delta_n+n^{-1}+\sqrt{\epsilon_n},
\]

where `epsilon_n` bounds the normalised initial optimisation error. The prediction ridge gradient contributes `O(n^{-1})` on the considered ball. On the event in (1), the population gradient in the selected span has norm at most `C nu_n`.

The weighted residual of any eligible column after projection on the selected span has norm bounded above and away from zero by the sparse design eigenvalue condition. Resolve `beta_star-beta_j` into this selected span and the residuals of its at most `s` remaining target columns. Convexity and population curvature imply

\[
 \max_{l\notin S_j}|\nabla R_0(\beta_j)^T b_{jl}|
       \geq c\sqrt{r_j/s}-C\nu_n, \tag{2}
\]

where `b_jl` is the corresponding normalised residual direction. To obtain (2), use `r_j<=grad R_0(beta_j)^T(beta_j-beta_star)`, bound the selected part by `C nu_n||beta_j-beta_star||_2`, and bound the sum of the remaining coefficients by `C sqrt(s)||beta_j-beta_star||_2`. Restricted curvature bounds that norm by `C sqrt(r_j)`.

The computed dual ordering attains the largest population residual gradient up to `C nu_n`. Bounded curvature along that direction supplies a population risk decrease at least its squared gradient divided by a fixed constant. Reoptimisation on the enlarged support has population excess over the population restricted optimum at most `C nu_n^2`: the approximate empirical optimality inequality and (1) give `kappa||Delta||_2^2/2<=C nu_n||Delta||_2+epsilon_n`. Thus

\[
 r_{j+1}\leq(1-c/s)r_j+C\nu_n^2,
 \qquad r_M\leq e^{-cM/s}r_0+Cs\nu_n^2. \tag{3}
\]

The paths stay in the required ball under the stated curvature condition. A direct localisation bound is available: nuisance stationarity bounds every fitted cluster effect by a constant; the profile loss then exceeds `c P_n m^{-1}sum_j|X_j^T beta|-P_n m^{-1}sum_j|Y_j|-C`. Uniform sparse small ball concentration bounds its first term below by `c'||beta||_2`. Every path fit has empirical objective no larger than the null fit up to its certified gap, so its coefficient norm is bounded in probability. Curvature is required on a ball containing this bound.

For any trained path candidate, its heldout loss difference from the target is bounded by `m^{-1}sum_j|X_j^T Delta|`. Conditional on training this has sub Gaussian scale `C||Delta||_2`. Population curvature then gives, simultaneously over the finite inner paths,

\[
 |(P_V-P)(\ell_0(\beta)-\ell_0(\beta^\star))|
 \leq a\{R_0(\beta)-R_0(\beta^\star)\}
       +C_a\log(CM/\delta)/n_V. \tag{4}
\]

The unobserved reference loss cancels when selecting the model size. Applying (4) to the selected size and the endpoint in (3) bounds the average selected population risk by `C[e^(-cM/s)+s nu_n^2+log(CM)/n]`. Each fold has a positive weight, so every selected inner coefficient has the corresponding bound.

The union of inner supports contains these fitted comparison vectors. Adding graph supports preserves that fact. The final outer training refit is an approximate empirical minimiser on the union. Apply (1) and the same restricted optimality inequality to its population restricted minimiser, whose risk is no larger than a comparison vector's risk. Hence, uniformly over the fixed requested targets and outer folds,

\[
 \|\widehat\beta_{fk}-\beta^\star\|_2=O_p(e_n),\quad
 e_n=e^{-cM/(2s)}+\sqrt{s}\nu_n+
                \sqrt{\log(CM)/n}. \tag{5}
\]

No equality between an inner selected path and a separate full training path is needed: the actual refitted support contains the comparison vectors.

## Positive curvature and learned contrasts

The cluster Schur variational identity gives

\[
 v^TH_bv=\min_a\{\|W_b^{1/2}(Xv-Za)\|_2^2+a^T\Lambda a\}\geq0.
\]

The project curvature moment proof applies to arbitrary unit design contrasts under the stated conditional sub Gaussian condition. In particular,

\[
 E|v^TH_bw|^r\leq r!C^r b^{1-r}\|v\|_2^r\|w\|_2^r,
 \qquad r\geq2. \tag{6}
\]

A coefficient net on each support and fixed mesh nets on the two unit spheres have total logarithmic size `O(DL)`. The deterministic profile Lipschitz bounds used in `thm:uniform-curvature` extend the coefficient net, and the usual sphere net inequality extends the bilinear bound. Bernstein's inequality applied to (6) yields the restricted operator bound

\[
 \sup_{\beta,A}\|(P_T-P)H_b(\beta)_{AA}\|_{op}
       =O_p(a_n),\quad
 a_n=\sqrt{DL/(nh)}+DL/(nh),\quad b\in\{h,h/2\}. \tag{7}
\]

Supports and coefficient radii in (7) include every considered selected set and its union with the sparse population contrast. Thus empirical restricted curvature is positive with probability tending to one when `a_n->0`.

Consider the standardised nodewise Lasso with sparse comparison contrast `c_0`. Its empirical gradient at that comparison is bounded by `C(a_n+e_n+h)` using (7), the population bandwidth expansion and the population curvature Lipschitz bound. At a penalty at least a sufficiently large multiple of this quantity, the Lasso basic inequality gives squared prediction error `C s_c(a_n+e_n+h)^2`. Its support size is at most `C s_c`. The support assertion follows from its KKT equations: on any subset of active coordinates the squared penalty is bounded by the sparse upper eigenvalue times prediction error. Taking a subset larger than a sufficiently large multiple of `s_c` contradicts this bound.

The ratio grid contains a suitable common candidate across the fixed inner folds and both bandwidths. Their standardised penalty anchors differ from the same population anchor by `O_p(a_n+e_n+h)`. If that population anchor is at most this order, the all zero graph has the asserted error. Otherwise the anchors are comparable, the geometric spacing is bounded, and the lower endpoint `n_training^(-1/2)` is below the required penalty order. A common ratio therefore gives the bound on all inner fits. The graph budget admits its supports by assumption.

Post Lasso quadratic refitting on the graph and outcome support retains the comparison prediction bound. For a heldout cluster, the difference of two quadratic nodewise losses has variance bounded by `C/h` times its population excess risk, by (6) and positive curvature. Bernstein's inequality and a union over the `O(M log n)` possible trained candidates give the validation oracle remainder `C log(CM log n)/(nh)`. The union over initial path sizes matters here: it permits their data selected size without treating it as independent of the inner validation data.

Each selected inner graph therefore contains a contrast with population L2 error bounded by

\[
 C\left[\sqrt{s_c}(a_n+e_n+h)+a_n+
           \sqrt{\log(CM\log n)/(nh)}\right].
\]

The outer graph union contains these comparison contrasts. Its two positive quadratic refits satisfy the same restricted optimality inequality using (7), with ridge contribution `O(n^{-1})`. Taking `2c_half-c_full` consequently gives

\[
 \|\widehat c_{fk}^R-c_0\|_2=O_p(r_n),\quad
 r_n=(1+\sqrt{s_c})(a_n+e_n+h+n^{-1})
          +\sqrt{\log(CM\log n)/(nh)}. \tag{8}
\]

For a fixed selected support, the population inverse block expansion additionally gives `2c_half-c_full=c_0,eta+O(h^2)` in Euclidean norm. Both inverses use the same reference ridge. Bound (8) also accounts for support approximation and is valid without asserting exact graph recovery.

## Common path expansion

Fix the fitted training coefficients and contrasts and define their population path criterion by expectation over a new independent cluster. At `theta_star=beta_star_k`, its coefficient error is

\[
 d_f=\widehat\beta_{fk}-\beta^\star
       -\widehat c_{fk}^R(\widehat\beta_{fk,k}-\beta_k^\star),
 \qquad d_{f,k}=0,\quad\|d_f\|_2=O_p(e_n).
\]

Let `M_b(t)` be the fold weighted population path score and `theta_b^P` its local zero. Since `c_0^T H_0 d_f=0`, a Taylor expansion gives

\[
 M_0(\theta^\star)=O_p(e_n^2+r_ne_n),\qquad
 \theta_0^P-\theta^\star=O_p(e_n^2+r_ne_n). \tag{9}
\]

The derivative of `M_0` tends to `J_0>0`. The density section expansion gives `M_b(t)=M_0(t)+b B(t)+O_p(b^2)` with the same fitted path throughout. Expanding at the zero of `M_0` and using the implicit function theorem yields

\[
 \theta_b^P=\theta_0^P-bB(\theta_0^P)/M_0'(\theta_0^P)+O_p(b^2).
\]

Consequently

\[
 2\theta_{h/2}^P-\theta_h^P-\theta^\star
           =O_p(h^2+e_n^2+r_ne_n). \tag{10}
\]

This expansion applies to the roots used in the implementation. Its single bandwidth directional derivatives remain positive.

## Inference theorem

Suppose the structural conditions hold and the errors in (5) and (8) tend to zero. Require

\[
 \sqrt n(h^2+e_n^2+r_ne_n)\longrightarrow0,\qquad
 \sqrt{DL}(e_n/\sqrt h+r_n)\longrightarrow0, \tag{11}
\]

and `DL sqrt(log n)/sqrt(n)->0`, `sqrt(D)r_n sqrt(L)->0`. Normalised numerical errors satisfy the displayed rules and are negligible in these rates. Then the paired root estimator and its recorded cluster standard error obey

\[
 \sqrt n(\widetilde\theta_k-\beta_k^\star)
 =-\frac1{J_0\sqrt n}\sum_{i=1}^n c_0^Tg_{i,0}(\beta^\star)+o_p(1), \tag{12}
\]

\[
 n\widehat{se}_k^2\ \longrightarrow_p\
 \sigma_k^2=\operatorname{Var}\{c_0^Tg_{i,0}(\beta^\star)\}/J_0^2,
 \quad(\widetilde\theta_k-\beta_k^\star)/\widehat{se}_k
       \ \Longrightarrow\ N(0,1). \tag{13}
\]

**Proof.** Conditional on the complementary training observations of one outer fold, its fitted coefficient and contrast are independent of the validation observations. The directional moment bound (6), integrated along coefficient segments, bounds the squared score substitution error by `C(e_n^2/h+r_n^2)`. When comparing the two bandwidth roots, use the deterministic population path through `beta_star` with direction `c_0` as reference at each bandwidth. Its displacement is `O(h)`. The extra displacement caused by contrast estimation is `O(h r_n)` and is included in the preceding bound. This choice of reference separates the bandwidth displacement from nuisance fitting error.

A sparse coefficient net of polynomial mesh and the deterministic profile Lipschitz bound give covering logarithm `O(DL)` for the local fitted score class. Bernstein's inequality bounds its centred empirical contribution, after multiplication by `sqrt(n)`, by

\[
 O_p\{\sqrt{DL}(e_n/\sqrt h+r_n)
                +DL\sqrt{\log n}/\sqrt n\}=o_p(1).
\]

To obtain the stated envelope term, centre the directions at `c_0`. The observed envelope is at most `C max_i max_j|X_ij^T c_0|+C sqrt(D)r_n max_ijk|X_ijk|`, of order `sqrt(log n)` under (11) and the last displayed rate condition. The polynomial net extension is negligible by choosing its mesh finer than the fixed powers in the profile Lipschitz bound. A union over the fixed outer folds is valid even though their training sets overlap. Independence is used separately for the validation observations of each fold.

The corresponding empirical directional derivatives converge to `J_0` by (7), population continuity and (5), (8). Taylor expansion of each monotone sample score about `theta_b^P` now gives

\[
 \sqrt n(\widehat\theta_b-\theta_b^P)
 =-J_0^{-1}n^{-1/2}\sum_i c_0^Tg_{i,0}(\beta^\star)+o_p(1)
\]

simultaneously for the two bandwidths. The unsmoothed reference is justified by the density section limit and uniform integrability of bounded dual scores times sub Gaussian design contrasts. Combining the two expansions with weights `(2,-1)` and using (10) proves (12). The independent cluster Lindeberg theorem applies because the influence variance is nondegenerate and its moments exceed order two.

For the variance, the fitted cluster contribution `V_i` differs from `-c_0^T g_i,0/J_0` in average squared norm by `o_p(1)`, using the same conditional score substitution bound, convergence of the two derivatives, and `e_n^2/h+r_n^2 -> 0`. Cauchy Schwarz controls the cross term in second moments. The iid law of large numbers and centring prove (13). The usual normal interval consequently has limiting coverage `1-alpha`. \(\square\)

## Explicit dimensional regime

For bounded `s`, `s_c` and `d_0`, take the implemented `M=ceil(log(n)^3)` asymptotically, `h=n^(-3/10)`, and numerical gaps at most order `n^(-3)` for initial paths and `n^(-4)` for final refits. If

\[
 \log(2np)=o\{n^{7/20}/(\log n)^3\}, \tag{14}
\]

then (1) has leading order `sqrt(DL/n)`, (5) has this same order, and (8) has order `sqrt(DL/(nh))+h`. The largest combined condition is `DL/sqrt(nh)->0`, which is exactly (14). Also `sqrt(n)h^2=n^(-1/10)->0`, while `h sqrt(DL)->0` and the envelope condition follow from (14). Thus (12) and (13) hold in a regime allowing `p` to grow exponentially in a positive power of `n`.

For growing sparsities the unsimplified rates (5), (8), (11) apply. They expose the cost of larger selected models and the number of required greedy steps. The inference theorem also applies to other training estimators satisfying its error conditions, including approximate sparse targets with their approximation errors added to `e_n` and `r_n`.

## Source connections

The entropy argument uses the polynomial sign counting result of Goldberg and Jerrum, *Machine Learning* 18, 131–148, 1995, Theorem 2.1 and Corollary 2.1. The remaining empirical process steps are the usual symmetrisation, covering and Bernstein arguments; the unbounded envelope can equivalently be treated by Adamczak, *Electronic Journal of Probability* 13, 1000–1034, 2008. The profile density expansion and curvature moments are derived in the project source cited above. Cross fitting with orthogonal scores is the established framework of Chernozhukov and colleagues. Here the common parameter path, paired root expansion, positive directional curvature and inner support union connect that framework to the implemented cluster quantile profile estimator.
