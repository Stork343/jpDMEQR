# Orthogonal Search for Clustered Quantile Profile Inference

Method identifier: `orthogonal_profile_score_20261025`.

The estimator targets the minimiser of the unsmoothed complete cluster profile risk with fixed positive definite working matrix `Lambda`. The implementation supports prespecified inference coordinates throughout the full design, including coordinates absent from the selected initial model.

## Prediction ridge and orthogonal coordinates

On a training sample of `n_T` independent clusters, set `w_ij=1/(n_T*m_i)` and let `A=diag(sqrt(w))*X`. For a selected set `S`, write `A_S=Q_S R_S` with orthonormal columns of `Q_S`. Define `B_S=diag(w^{-1/2})Q_S` and `a=R_S beta_S`. Jointly minimise

\[
 \frac1{n_T}\sum_i\left\{\frac1{m_i}\sum_j
 \rho_\tau(Y_{ij}-B_{ij,S}^Ta-Z_{ij}^T\gamma_i)
 +\frac12\gamma_i^T\Lambda\gamma_i\right\}
 +\frac{\|a\|_2^2}{2n_T}.
\]

The small ridge equals `beta_S^T Sigma_hat_S beta_S/(2*n_T)` with `Sigma_hat=X^T diag(w) X`. It therefore regularises the fitted linear predictor and is invariant to the chosen orthonormal basis for the selected span. Its coefficient tends to zero with the number of training clusters. The inferential profile score contains neither this initial ridge nor a variable selection penalty.

The primal dual solver operates on the selected coordinates and all cluster nuisance vectors simultaneously. The conjugate of `||a||^2/2` is `||v||^2/2`, so every box feasible dual vector supplies a computable gap for each selected convex problem.

## Conditional variable contribution

Let `Q` be the current weighted orthogonal basis. For unselected column `j`, put

\[
 v_j=(I-QQ^T)A_j,\qquad d_j=\|v_j\|_2,\qquad
 b_j=\operatorname{diag}(w^{-1/2})v_j/d_j.
\]

Select the coordinate maximising `abs(b_j^T u)`, where `u` is the current residual dual vector. This measures the candidate contribution after removing variation already represented by the selected model. A zero residual design norm denotes a redundant extension.

The basis gains the column `v_j/d_j`, and the triangular factor is updated. The preceding basis columns remain unchanged. The enlarged selected problem uses a zero initial coefficient on the new basis column and reuses the nuisance and dual iterates. Modified Gram Schmidt with reorthogonalisation maintains the weighted design representation. Original coefficient units are recovered by triangular solution.

## Extension bound

For the current basis, the dual objective is

\[
 D_S(u)=Y^Tu-\tfrac12\|B_S^Tu\|_2^2
 -\tfrac12\sum_i(Z_i^Tu_i)^T\Lambda^{-1}(Z_i^Tu_i).
\]

Adding the orthogonal candidate `b_j` changes this feasible dual value by exactly `-(b_j^T u)^2/2`. If the current normalised primal dual gap is `epsilon_S`, the attainable decrease in the normalised criterion is at most

\[
 \epsilon_S+\frac{(b_j^Tu)^2}{2n_T}.
\]

The selection rule orders these bounds. It uses the same coordinates in the search statistic and in the subsequent optimisation. The observed objective decrease and the bound are independently checked in `validate_orthogonal.R`.

## Model size and cross fitting

Five whole cluster folds build complete search paths. Always included covariates are supplied by zero inclusion indicators; one denotes an eligible search variable. The candidate path size is bounded by the smaller training fold's `floor(n_T/log(n_T))`, after allowing for the always included coordinates and the available columns.

The selected size minimises the total exact heldout profile loss. All training transformations and search decisions use their corresponding training subjects. Heldout profile loss optimises a nuisance vector for each validation cluster, which matches the population profile target. The selected heldout coefficient vectors are reused for inference.

Each cluster's paired score and positive single bandwidth curvature are evaluated at the coefficient fitted without its fold. All `p` original columns enter these calculations. The positive curvature factorisation, requested nodewise directions and projected score inversion are defined in GREEDY_METHOD.md. The nodewise path endpoint is `min(0.05,n^(-1/2))` relative to its observed path anchor, allowing its smallest penalty to decrease with sample size.

## Low dimensional nuisance computation

The unsmoothed random intercept problem is solved by a sorted scalar stationarity sweep. For an intercept and slope, whitening the ridge and restricting to residual zero lines gives a finite collection of scalar problems and adjacent cell candidates. Every accepted point satisfies a complete primal dual gap check. FINITE_PROFILE.md gives the completeness argument and describes degenerate ties.

The smoothed nuisance evaluations used for score inference have their own gradient tolerance and exact Schur identity. Each cluster requires a problem with nuisance dimension `q`.

## Statistical results and remaining derivations

The validation risk inequality in GREEDY_METHOD.md applies to these training paths because each candidate is determined by its training clusters. Heldout risk uses the complete cluster profile criterion. Population curvature converts the resulting path approximation bound into coefficient error bounds.

The conditional score inference argument in THEORY.md and GREEDY_METHOD.md requires initial error and inverse direction rates, the profile bandwidth expansion and a local identification condition. The algorithm supplies numerical gaps, observed direction defects and local score derivatives. Establishing useful approximation rates for the actual orthogonal search path and data selected direction parameters remains a mathematical research requirement. Signal and design conditions enter that primitive rate analysis.

## Complexity

The weighted design and its projections are retained. Each added coordinate requires a complete design score product and a projection update, costing order `N*p`. Each optimisation iteration then operates on `d` selected coordinates, costing order `N*d+N*q+n*q^2`. The path cost is

\[
 O\!\left(MNp+\sum_{d=1}^M I_d(Nd+Nq+nq^2)\right),
\]

where `I_d` is the measured iteration count. Selected fits retain duality gap certificates. Curve factors and nodewise regressions use rectangular designs, and score inversion caches the high dimensional projections before scalar root solving.

## Interpretation of innovation

Orthogonal greedy search, primal dual optimisation, nodewise regression and cross fitting have published precedents. The research-specific constructions are the finite complete quantile profile solver, the prediction ridge extension bound in the evolving orthogonal basis, and the combination with positive profile curvature and paired cluster score inversion. Their numerical implementation and algebraic identities have direct checks. Priority and a complete primitive statistical theorem require comparison with the closest literature and further derivation.
