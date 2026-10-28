source("research/submission_completion_20261027/load_comparators.R")
records <- list()
for (p in c(30L, 80L, 250L)) for (tau in c(.25, .5, .75)) {
  set.seed(202610280L + p + as.integer(tau * 100))
  X <- cbind("(Intercept)" = 1, matrix(rnorm(300L * (p - 1L)), 300L))
  colnames(X)[-1L] <- paste0("x", seq_len(p - 1L))
  y <- 1 + X[, 2L] - .5 * X[, 3L] + rnorm(300L)
  penalty <- sqrt(tau * (1 - tau) * 300L * log(p)) * c(0, rep(1, p - 1L))
  started <- proc.time()[[3L]]
  reference <- fit_quantile_l1_block_v2(X, y, tau, penalty)
  reference_seconds <- proc.time()[[3L]] - started
  started <- proc.time()[[3L]]
  candidate <- submission_quantile_column_fit(X, y, tau, penalty)
  candidate_seconds <- proc.time()[[3L]] - started
  difference <- max(abs(reference$coefficients - candidate$coefficients))
  reference_objective <- sum(reference$residuals * (tau - (reference$residuals < 0))) +
    sum(penalty * abs(reference$coefficients))
  objective_difference <- abs(reference_objective - candidate$objective)
  if (difference > 2e-5 || objective_difference > 1e-6 * max(1, reference_objective))
    stop("The quantile regression solution changed after column generation.")
  records[[length(records) + 1L]] <- data.frame(p = p, tau = tau, coefficient_difference = difference,
    objective_difference = objective_difference, maximum_dual_violation = candidate$maximum_dual_violation,
    primal_dual_gap = candidate$primal_dual_gap, reference_seconds = reference_seconds, candidate_seconds = candidate_seconds)
  write.csv(do.call(rbind, records), "results/scratch/submission_completion_20261027/quantile_column_checks.csv", row.names = FALSE)
}
print(do.call(rbind, records))
