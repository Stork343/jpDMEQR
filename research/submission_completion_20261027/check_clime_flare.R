source("research/submission_completion_20261027/load_comparators.R")
source("research/submission_completion_20261027/clime_flare.R")
records <- list()
for (rank in c(20L, 40L)) {
  set.seed(202610270L + rank)
  X <- matrix(rnorm(rank * 30L), rank, 30L)
  H <- crossprod(X) / rank
  dimnames(H) <- list(paste0("x", 1:30), paste0("x", 1:30))
  for (gamma in c(.2, .4, .6)) {
    reference <- solve_clime_matrix_iid_v2(H, gamma)
    candidate <- solve_clime_matrix_iid_v2(H, gamma, solver_opts = list(algorithm = "flare_warm"))
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
  package <- flare::sugm(H, lambda = .4, method = "clime", standardize = FALSE, perturb = FALSE,
    verbose = FALSE, prec = 1e-5, max.ite = 10000L)
  optimized <- submission_flare_matrix(H, .4)
  difference <- max(abs(package$icov1[[1L]] - optimized$matrices[[1L]]))
  if (difference > 1e-7) stop("The flare matrix algebra changed its CLIME output.")
}
write.csv(do.call(rbind, records), "results/scratch/submission_completion_20261027/clime_flare_checks.csv", row.names = FALSE)
print(do.call(rbind, records))
