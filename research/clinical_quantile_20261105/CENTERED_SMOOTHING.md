# Centered smoothing for clinical quantile estimation

The scientific parameter and the unsmoothed patient orthogonal score are those in METHOD.md. The preliminary regression uses a centred smooth loss to reduce computation with dense protein coefficients.

Let H_delta(r) equal r squared divided by 2 delta for |r|<=delta and |r|-delta/2 otherwise. Write ell_(tau,delta)(r)=tau H_delta(r) for r>=0 and (1-tau)H_delta(r) for r<0. Set a_ij=N/(n m_i), delta=n^(-1/5) and s=(2 tau-1)delta/2. The preliminary estimator minimises

\[
N^{-1}\sum_{i,j}\ell_{\tau,\delta}
\{a_{ij}(Y_{ij}-W_{ij}^Tb)+s\}
+\lambda\{\alpha\|b_P\|_1+(1-\alpha)\|b_P\|_2^2/2\}.
\]

The derivative of ell is a smoothed check score whose kernel is a mixture of Uniform(-delta,0) with weight 1-tau and Uniform(0,delta) with weight tau. This kernel has mean s. Shifting the residual by s centres the kernel. Its variance is delta squared times {1/3-(2 tau-1)^2/4}. For a conditional distribution with a differentiable density, the expectation of the centred score at the true quantile is O(delta squared). Bounded cluster sizes keep the row scales a_ij bounded above and away from zero. Positive local curvature transfers this score error to the population coefficient approximation.

The smoothing contribution to a preliminary prediction rate is therefore O(n^(-2/5)). It satisfies sqrt(n) times its square ->0. Products with estimated density and instrument errors still have to satisfy the conditions in METHOD.md. The patient orthogonal update uses the original indicator score, and its target remains the stated conditional quantile coefficient.

The solver uses skglm QuantileHuber and WeightedL1_plus_L2 with AndersonCD. Clinical columns are transformed to an orthogonal basis, and protein columns are residualised against that basis. Write theta_C for the coefficients of that clinical basis. The objective additionally contains (||theta_C|| squared+||b_P|| squared)/(2n). This vanishing ridge defines a unique solution, including when a clinical indicator occurs only in observations in the linear part of the loss. All clinical directions remain free of L1 selection. The n^(-1) ridge is included in the numerical objective and KKT checks; its population contribution vanishes faster than n^(-1/2) for bounded comparison coefficients.

Returned coefficients are mapped back to the original training coordinates. The quadratic ridge in original coordinates is computed with the clinical projection matrix and the protein coefficient norm. Full original-coordinate KKT residuals are recomputed for every path entry.

The skglm source is pinned to commit c73abe71b89752ec9c6444932fffb5bb73ff7a9c. The weight capable penalty is present in that source version. The package is installed in the project environment `results/scratch/cq_solver_env`. Direct mathematical checks compare the smooth objective with an independently formulated quadratic programme on actual training data. Deterministic population integration checks the centred kernel and the approximation order.
