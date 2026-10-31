library(data.table)
source("research/score_validity_20261026/load_joint.R")
source("research/milk_application_20261031/weighted_solver.R")
cohort <- readRDS("results/raw/application/milk_scs_20261031/cohorts.rds")$eligible
data <- sb_data(cohort$y, cohort$X[, 1:13], cohort$Z, cohort$cluster_id)
initial <- drop(qr.solve(data$X * sqrt(data$weight), data$y * sqrt(data$weight)))
checks <- list()
for (tau in c(0.25, 0.5, 0.75)) {
  for (multiplier in c(1, 0.5)) {
    h <- data$n^(-0.3) * multiplier
    original <- milk_reference_fit_small(data, initial, tau, h, matrix(1, 1L, 1L))
    stable <- milk_fit_small(data, initial, tau, h, matrix(1, 1L, 1L))
    errors <- c(coefficient = max(abs(original$beta - stable$beta)),
      influence = max(abs(original$projected_scores - stable$projected_scores)),
      objective = abs(original$objective - stable$objective),
      hessian = max(abs(original$hessian - stable$hessian)),
      original_gradient = stable$kkt, prediction_identity = stable$prediction_error)
    tolerances <- c(2e-6, 2e-5, 1e-9, 2e-5, 1e-8, 1e-7)
    if (any(!is.finite(errors)) || any(errors > tolerances)) stop("The observed-data basis check failed.")
    checks[[length(checks) + 1L]] <- data.table(tau = tau, h = h, check = names(errors),
      error = as.numeric(errors), tolerance = tolerances)
  }
}
fwrite(rbindlist(checks), "results/preflight/milk_scs_20261031/weighted_solver_checks.csv")
print(rbindlist(checks)[, .(maximum_error = max(error)), by = check])
