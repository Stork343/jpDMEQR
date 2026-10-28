source("research/submission_completion_20261027/load.R")
records <- list()
for (rank in c(20L, 40L)) {
  set.seed(202610270L + rank)
  X <- matrix(rnorm(rank * 30L), rank, 30L)
  H <- crossprod(X) / rank
  dimnames(H) <- list(paste0("x", 1:30), paste0("x", 1:30))
  for (gamma in c(.2, .4, .6)) {
    reference <- solve_clime_matrix_iid_v2(H, gamma)
    candidate <- solve_clime_matrix_iid_v2(H, gamma, solver_opts = list(basis_policy = "reuse_primal_feasible"))
    if (!identical(reference$feasible, candidate$feasible)) stop("CLIME feasibility differs across simplex algorithms.")
    error <- if (reference$feasible) max(abs(rowSums(abs(reference$omega_raw)) - rowSums(abs(candidate$omega_raw)))) else 0
    residual <- if (candidate$feasible) max(abs(candidate$omega_raw %*% H - diag(30))) else NA_real_
    if (error > 2e-6 || (candidate$feasible && residual > gamma + 1e-6)) stop("CLIME objective or residual verification failed.")
    records[[length(records) + 1L]] <- data.frame(rank = rank, gamma = gamma, feasible = candidate$feasible,
      objective_error = error, maximum_residual = residual)
  }
}
write.csv(do.call(rbind, records), "results/scratch/submission_completion_20261027/clime_basis_checks.csv", row.names = FALSE)
print(do.call(rbind, records))
