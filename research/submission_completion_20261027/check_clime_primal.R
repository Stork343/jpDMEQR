source("research/submission_completion_20261027/load.R")
source("research/submission_completion_20261027/clime_primal.R")
records <- list()
for (rank in c(20L, 40L)) {
  set.seed(202610270L + rank)
  X <- matrix(rnorm(rank * 30L), rank, 30L)
  H <- crossprod(X) / rank
  dimnames(H) <- list(paste0("x", 1:30), paste0("x", 1:30))
  for (gamma in c(.2, .4, .6)) {
    reference <- solve_clime_matrix_iid_v2(H, gamma)
    candidate <- solve_clime_matrix_iid_v2(H, gamma, solver_opts = list(algorithm = "primal_warm"))
    if (!identical(reference$feasible, candidate$feasible)) stop("CLIME feasibility differs across formulations.")
    error <- if (reference$feasible) max(abs(reference$omega_raw - candidate$omega_raw)) else 0
    objective_error <- if (reference$feasible)
      max(abs(rowSums(abs(reference$omega_raw)) - rowSums(abs(candidate$omega_raw)))) else 0
    residual <- if (candidate$feasible) max(abs(candidate$omega_raw %*% H - diag(30))) else NA_real_
    if (objective_error > 2e-6 || (candidate$feasible && residual > gamma + 1e-6))
      stop("CLIME objective or residual verification failed.")
    records[[length(records) + 1L]] <- data.frame(rank = rank, gamma = gamma, feasible = candidate$feasible,
      coefficient_error = error, objective_error = objective_error, maximum_residual = residual)
  }
}
write.csv(do.call(rbind, records), "results/scratch/submission_completion_20261027/clime_primal_checks.csv", row.names = FALSE)
print(do.call(rbind, records))
