# Certified selected fits for paired observations

For a selected weighted orthogonal design $X$, the unsmoothed search problem is

\[
 \min_{\theta,\gamma}
 \frac12\|\theta\|^2+\frac{\lambda}{2}\|\gamma\|^2
 +\sum_{i,j}\frac12\rho_\tau(y_{ij}-x_{ij}^T\theta-\gamma_i).
\]

Dividing by the number of lines gives the registered complete profile objective and prediction ridge. All lines in this application have exactly two observations, and the nuisance design is an intercept.

Introduce a residual vector and apply two-block ADMM. The residual update is the exact proximal map of the check loss. For the quadratic update, let $W=I/2$ and let $C$ contain the within-line averages of the design rows as columns. Eliminating the line effects gives

\[
 H_\rho=I+\rho X^TWX-\frac{\rho^2}{\lambda+\rho}CC^T.
\]

This matrix has only the selected design dimension. Its Cholesky factor is reused throughout a selected fit. The implementation uses fixed augmentation parameter one and overrelaxation 1.6.

When the residual active pattern permits it, an equality-constrained quadratic calculation refines the solution. A line with one zero residual determines its intercept from that observation. Two zero residuals contribute a linear equality on the selected coefficients. The resulting system has at most twice the selected dimension. A proposed refinement is accepted only when it meets the complete primal-dual certificate.

The dual scores satisfy $s_{ij}\in[(\tau-1)/2,\tau/2]$. Their objective is

\[
 D(s)=y^Ts-\frac12\|X^Ts\|^2
      -\frac1{2\lambda}\sum_i(s_{i1}+s_{i2})^2.
\]

Every returned fit recomputes the primal objective using its actual residuals and compares it with this box-feasible dual value. The normalised gap must satisfy the original requested tolerance. This certification covers both an ADMM iterate and an active-face refinement. Exhausting the iteration budget terminates the attempted fit.

The source model, variable search, budgets and validation criterion are unchanged. The implementation is verified against the independent Clarabel formulation on actual SoyNAM observations, with recorded comparisons in `results/preflight/soynam_20261031/pair_qp_checks.csv`.
