# SoyNAM pooled quantile Lasso

The prediction comparator minimises the complete quantile Lasso objective

\[
 N^{-1}\sum_{j=1}^N\rho_\tau(y_j-x_j^T\beta)
 +\lambda\sum_{k\in G}\sigma_k|\beta_k|,
\]

where $G$ contains all genomic columns and $\sigma_k$ is their training standard deviation. Prespecified canopy, family and field adjustment terms are unpenalised. The two observations of each line make equal-row and equal-line fitting weights identical. The predictor is evaluated on complete held-out lines.

Centred and scaled genomic coordinates preserve this objective, with the induced intercept represented by the exhaustive block indicators. In the transformed design $X_s$, the dual programme maximises $y^Ta$ subject to $a_j\in[\tau-1,\tau]$ and $|(X_s^Ta)_k|\le N\lambda$ for genomic columns. The unpenalised columns have equality constraints.

An optimal dual solution of the adjustment-only fit supplies the entry penalty $\lambda_{max}=\max_{k\in G}|(X_s^Ta)_k|/N$. Each of three internal training samples uses the same 40 geometric ratios from 1 to 0.01. Selection uses their complete-line validation pinball loss, and the final training model uses the selected ratio.

The compiled solver calls HiGHS 1.15.1 through its official C++ API. It retains one LP object and its simplex basis through successive penalties. At every solve, it evaluates all genomic dual constraints. Violating constraints are added in descending violation order. Redundant basic inequality rows may be removed between penalties; complete pricing is repeated at the next point. Acceptance requires the original primal and dual objectives to agree within $10^{-7}$ per observation, the full dual violation to be at most $10^{-8}$, and the coefficient transformation to preserve predictions within $10^{-9}$.

`check_pooled_compiled.R` verifies three quantiles and three penalty levels, including zero, against independent `quantreg::rq.fit.fnb` solutions on actual field observations. A complete 40-point path with all 8520 prediction columns is also checked against an independently constructed sequence of R/HiGHS programmes. Every internal and final training path retains all coefficients and numerical certificates. This comparator provides predictions and loss comparisons.

The quantile Lasso formulation follows the standard penalised check-loss construction used in the repository's POOLED-QR-LASSO benchmark. The source reference is Belloni and Chernozhukov (2011), DOI 10.1214/10-AOS827. The numerical certificate follows directly from the displayed primal and dual programmes.
