# Paired population computation

## Full design streamed calculation

`R/v2_simulation_v2_part05.R` evaluates the same complete ridge-profile criterion as `profile_bandwidth_pair_v2`. It retains every fixed-effect column. Each nuisance problem is solved to an explicitly recorded gradient tolerance at both bandwidths, using the same observations and Lambda.

For one bandwidth, let `W_i=diag(phi_i)/m_i`, `U_i=Z_i'W_i X_i`, and `B_i=Z_i'W_i Z_i+Lambda`. If `R_i=chol(B_i)` is upper triangular, the exact Hessian sum is

\[
\sum_i X_i^{\mathsf T}W_iX_i-\sum_i T_i^{\mathsf T}T_i,
\qquad T_i=(R_i^{\mathsf T})^{-1}U_i.
\]

The q=1/2 formulas evaluate this factorisation directly. The second identity uses `X_tilde_i'W_i X_tilde_i+A_i'Lambda A_i`. Both expressions are accumulated in batches and compared. The batch calculation constructs no per-cluster p by p matrices.

Cluster score sums and uncentred crossproducts are accumulated for the half and full bandwidths. Their means are subtracted after division by the total number of clusters. The paired score crossproduct is also computed directly and compared with the two component crossproducts and their cross-bandwidth product. Each target optimiser evaluates score and objective without constructing covariance matrices; full geometry and moments are computed at the accepted root.

`src/population_products.cpp` provides the registered two-argument `.Call` routine `jpDMEQR_population_crossprod`. The explicit backend `Eigen_single_thread` uses Eigen matrix products. `R_BLAS` remains available for reference comparisons. Each backend uses its own floating point summation order, and the numerical differences are checked on identical data.

## Conditional Gaussian design integration

This calculation is used only when `profile_pair_gaussian_row_conditions_v2` accepts the frozen scenario. At the generated structural coefficient, the residual before nuisance fitting depends on cluster size, nuisance variables and errors. The fixed-effect rows are independent centred Gaussian vectors with covariance Sigma_x, independent of those residual variables and independent across observations.

The eligible registered direction cases are P01, P05, P06, A01 through A12, and D01 through D06. The P and A cases have Gaussian AR1 fixed-effect rows independent of the nuisance variables. A10 through A12 change the support size, which preserves this independence. D01 uses independent coordinates, D02 changes the AR1 correlation, D03 uses Gaussian block factors, D04 applies a deterministic precision transformation to Gaussian rows, D05 combines independent Gaussian factors and errors within each row, and D06 uses a Gaussian common component within each row. Each D design has independent rows and the explicitly computed full covariance Sigma_x.

Write `L_i=W_i-W_i Z_i B_i^-1 Z_i'W_i`. Conditional on the nuisance variables,

\[
E_X(X_i^{\mathsf T}L_iX_i)=\operatorname{tr}(L_i)\Sigma_x,
\qquad
\operatorname{tr}(L_i)=\operatorname{tr}(W_i)
-\operatorname{tr}\{B_i^{-1}Z_i^{\mathsf T}W_i^2Z_i\}.
\]

The conditional mean of the paired cluster score is exactly zero. Its covariance is

\[
\frac1{m_i^2}\sum_j
\{2\psi_{h/2}(r_{ij,h/2})-\psi_h(r_{ij,h})\}^2\Sigma_x.
\]

The half/full cross covariance replaces the square by the corresponding product of the two component scores. The paired Hessian uses twice the half-bandwidth trace minus the full-bandwidth trace. All averages use the requested 100000 independent clusters, and four independently seeded repetitions provide the between-repeat MCSE. The full p by p matrices and all requested directions are retained.

The generator consumes the original design-normal draws in their original order. This preserves the reference draws for cluster sizes, random effects, visit times and errors under the same seed, while Gaussian X is integrated conditionally. The comparison scripts save the difference between finite full-X Monte Carlo matrices and conditional-integral estimates, together with its MCSE. Those finite matrices have different Monte Carlo variation.

The conditional engine rejects correlated X and random effects, informative-size settings, nonlinear omitted nuisance terms, heteroskedastic errors, non-Gaussian errors, and mismatched working nuisance designs. C06, C07, C14, C15 and C16 therefore use full-design streamed target optimisation and independent full-design streamed direction calculations.

## Execution identity and outputs

The construction version is `profile-pair-population-assets-20261022-v3`. Assets are stored under `results/intermediate/profile_pair_20261022/construction_v3/targets` and `directions`. The source and evaluation point are explicit, including the structural-definition reference and the independent population target reference. Every repeat retains nuisance, profile identity, curvature and covariance diagnostics.

The formal entry is `scripts/research/build_profile_pair_population_final_20261022.R`. It requires a paired versioned registry, the current implementation commit, a new run identifier and an installed package containing the registered native routine. It fixes `n_population=100000` and `repeats=4`. Default concurrency is two cases, with one thread per matrix calculation. Original assets and run directories remain available.

An asset retains its actual construction commit. A subsequent manuscript or comparator commit may read it only when that construction commit is an ancestor and every population dependency listed by `profile_pair_population_source_paths_v2()` is unchanged between the two commits. The generating design, evaluation point, bandwidth and construction version are still validated separately. Changing any population dependency requires newly constructed assets.
