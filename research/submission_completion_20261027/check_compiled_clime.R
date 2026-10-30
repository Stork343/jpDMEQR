source("research/submission_completion_20261027/load_comparators.R")
source("research/submission_completion_20261027/load_compiled_clime.R")
records <- list()
for (rank in c(20L, 40L)) {
  set.seed(202610270L + rank)
  X <- matrix(rnorm(rank * 30L), rank, 30L)
  H <- crossprod(X) / rank
  dimnames(H) <- list(paste0("x", 1:30), paste0("x", 1:30))
  for (gamma in c(.2, .4, .6)) {
    reference <- solve_clime_matrix_iid_v2(H, gamma)
    candidate <- solve_clime_matrix_iid_v2(H, gamma, solver_opts = list(algorithm = "compiled_certified"))
    if (!identical(reference$feasible, candidate$feasible)) stop("CLIME feasibility differs between solvers.")
    error <- objective_error <- residual <- dual_residual <- gap <- loss_error <- NA_real_
    if (candidate$feasible) {
      error <- max(abs(reference$omega_raw - candidate$omega_raw))
      objective_error <- max(abs(rowSums(abs(reference$omega_raw)) - rowSums(abs(candidate$omega_raw))))
      dual <- do.call(rbind, lapply(candidate$raw_grid_records, function(x) x$solver_records[[1L]]$dual))
      residual <- max(abs(candidate$omega_raw %*% H - diag(30)))
      dual_residual <- max(abs(dual %*% H))
      gap <- max(abs(rowSums(abs(candidate$omega_raw)) - diag(dual) + gamma * rowSums(abs(dual))))
      ordinary_loss <- sum((H %*% candidate$omega_symmetric - diag(30))^2) / 30
      loss_error <- abs(ordinary_loss - submission_sqr_sparse_loss_cpp(H, candidate$omega_symmetric))
      complete <- submission_clime_highs_cpp(H, gamma, FALSE, H, Inf)
      bounded <- submission_clime_highs_cpp(H, gamma, FALSE, H, ordinary_loss / 2)
      reference_columns <- colSums((H %*% candidate$omega_symmetric - diag(30))^2) / 30
      evaluated <- which(!is.na(bounded$column_losses))
      if (complete$pruned || complete$evaluated_columns != 30L ||
          abs(complete$loss_lower_bound - ordinary_loss) > 1e-10 ||
          !bounded$pruned || bounded$loss_lower_bound <= ordinary_loss / 2 ||
          bounded$loss_lower_bound > ordinary_loss + 1e-10 ||
          (length(evaluated) && max(abs(bounded$column_losses[evaluated] - reference_columns[evaluated])) > 1e-10) ||
          any(bounded$column_lower_bounds > reference_columns + 1e-10))
        stop("The partial-column loss did not certify the complete validation comparison.")
      if (error > 2e-5 || objective_error > 2e-6 || residual > gamma + 1e-7 ||
          dual_residual > 1 + 1e-7 || gap > 2e-6 || loss_error > 1e-10)
        stop("The independent CLIME certificate or validation-loss comparison failed.")
    }
    records[[length(records) + 1L]] <- data.frame(rank = rank, gamma = gamma, feasible = candidate$feasible,
      coefficient_error = error, objective_error = objective_error, primal_residual = residual,
      dual_residual = dual_residual, dual_gap = gap, validation_loss_error = loss_error)
  }
}
output <- "results/preflight/sqr_acceleration_20261030"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
write.csv(do.call(rbind, records), file.path(output, "compiled_clime_checks.csv"), row.names = FALSE)
print(do.call(rbind, records))
box_records <- list()
for (rank in c(4L, 5L)) for (radius in c(.05, .4)) {
  set.seed(202610300L + rank)
  C <- matrix(rnorm(12L * 5L), 12L, 5L)
  if (rank == 4L) C[, 5L] <- C[, 1L]
  residual <- rnorm(12L)
  radii <- rep(radius, 5L)
  lower <- submission_sqr_box_bound_cpp(residual, C, radii)
  model <- osqp::osqp(P = Matrix::Matrix(2 * crossprod(C), sparse = TRUE),
    q = drop(2 * crossprod(C, residual)), A = Matrix::Diagonal(5L), l = -radii, u = radii,
    pars = osqp::osqpSettings(verbose = FALSE, eps_abs = 1e-9, eps_rel = 1e-9, max_iter = 10000L))
  fitted <- model$Solve()
  if (fitted$info$status != "solved") stop("The independent box-programme comparison did not converge.")
  minimum <- sum((residual + C %*% fitted$x)^2)
  if (!is.finite(lower) || lower <= 0 || lower > minimum + 1e-8)
    stop("The validation certificate exceeds the independent box-programme optimum.")
  box_records[[length(box_records) + 1L]] <- data.frame(rank = rank, radius = radius,
    certified_lower_bound = lower, independent_minimum = minimum)
}
write.csv(do.call(rbind, box_records), file.path(output, "box_loss_checks.csv"), row.names = FALSE)
print(do.call(rbind, box_records))
