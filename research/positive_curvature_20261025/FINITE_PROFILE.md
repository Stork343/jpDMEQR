# Finite Profile Evaluation with One or Two Nuisance Coordinates

For a cluster of size `m`, the unsmoothed problem is

\[
 \min_\gamma\;m^{-1}\sum_j\rho_\tau(y_j-z_j^T\gamma)
                         +\tfrac12\gamma^T\Lambda\gamma.
\]

Its dual maximises

\[
 y^Tu-\tfrac12(Z^Tu)^T\Lambda^{-1}(Z^Tu),
 \qquad (\tau-1)/m\le u_j\le\tau/m.
\]

Every feasible dual vector supplies a lower bound. The implementation accepts a candidate only when the complete primal dual gap is at most `1e-10*max(1,abs(primal))`. Box feasibility is checked before computing that certificate.

## One coordinate

For nonzero `z_j`, put `t_j=y_j/z_j`, `w_j=abs(z_j)/m`, and let `tau_j=tau` for positive `z_j` and `tau_j=1-tau` for negative `z_j`. Away from the sorted knots, the derivative is

\[
 \lambda\gamma+\sum_{t_j<\gamma}w_j-\sum_jw_j\tau_j.
\]

Its root on an interval is explicit. At a knot, the left and right derivatives determine whether that knot is the minimiser. Sorting and a single sweep find the unique optimum. Zero design entries contribute a constant loss. For zero residuals, the algorithm assigns dual contributions within their box intervals to satisfy `lambda*gamma=sum(z_j*u_j)`. Tied knots and signed design entries are covered by the same calculation.

The cost is order `m*log(m)`, followed by a sweep over at most `m` distinct knots.

## Two coordinates

Write `Lambda=R^T R`, `theta=R gamma`, and `A=Z R^{-1}`. The criterion becomes

\[
 f(\theta)=m^{-1}\sum_j\rho_\tau(y_j-a_j^T\theta)+\tfrac12\|\theta\|_2^2.
\]

Residual zero lines partition the plane into cells. In a cell the loss is quadratic with Hessian equal to the identity. Its stationary point is `theta=A^T u`, where each nonzero residual fixes `u_j` at the corresponding box endpoint.

On line `a_j^T theta=y_j`, write `theta=o_j+v_j*t`, where `o_j` is the perpendicular projection of zero onto the line and `v_j` is a unit tangent. Since `o_j^T v_j=0`, the restricted quadratic is `t^2/2` plus a constant. The restricted loss is a one dimensional ridge quantile problem with responses `y-A o_j` and design `A v_j`. The preceding scalar sweep therefore computes its unique minimiser.

At a restricted minimiser, nonzero residuals fix their dual endpoints. The remaining box variables satisfy two stationarity equations. A feasible polytope with two independent equalities has a vertex with at most two free coordinates. Enumerating those free pairs and endpoint assignments finds a feasible dual vector whenever one exists. Collinear active rows reduce to scalar interval allocation. The algorithm also examines the adjacent endpoint patterns and their cell stationary points.

## Completeness

If the global minimiser lies on a residual zero line, it must minimise the convex restriction to that line. The line search therefore includes it, and dual feasibility holds by the optimality conditions.

If the global minimiser is in a cell interior, it equals that cell's stationary point. Choose a residual zero line closest to it in the whitened Euclidean metric. The perpendicular projection onto this line belongs to the closure of the same cell because the open ball of that radius contains no residual zero line. At the projection, the tangential derivative of the cell quadratic is zero. Any kink on the line changes only its normal subgradient, so the projection minimises the whole convex line restriction. One of the adjacent endpoint patterns is the pattern of the cell containing the global minimiser. Its stationary point is therefore enumerated.

These two cases cover the unique minimiser. Coincident residual lines are handled by the active dual feasibility calculation and endpoint assignments. The primal dual certificate verifies the numerical candidate independently of this enumeration argument.

With generic outcomes, at most two residuals are active at a candidate and the work is order `m^2*log(m)`. Degenerate ties can require endpoint enumeration with a factor exponential in `m`. The intended regime has a fixed small maximum cluster size. The implementation retains the finite dimension cap on exceptional coincident active sets and reports an error if that cap is exceeded.

## Numerical verification

The validation script compares 300 cluster problems with an independently assembled CLARABEL quadratic programme. The cases include cluster sizes two through twelve, non diagonal working penalties, duplicate outcomes, a rank one nuisance design, zero design rows and nearly coincident visit times. Objective differences and coefficient differences are saved in the machine readable validation output.
