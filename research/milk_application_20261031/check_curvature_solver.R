library(data.table)
source("research/score_validity_20261026/load_joint.R")
source("research/milk_application_20261031/nodewise.R")
source("research/milk_application_20261031/weighted_solver.R")
cohort <- readRDS("results/raw/application/milk_scs_20261031/cohorts.rds")$eligible
columns <- unique(c(1:13, as.integer(seq(14, ncol(cohort$X), length.out = 30))))
data <- sb_data(cohort$y, cohort$X[, columns], cohort$Z, cohort$cluster_id)
initial <- numeric(ncol(data$X))
initial[2:13] <- median(data$y)
checks <- list()
for (tau in c(0.25, 0.5, 0.75)) {
  for (h in data$n^(-0.3) * c(1, 0.5)) {
    profile <- ds_profile(data, initial, tau, h, matrix(1, 1L, 1L))
    factor <- sb_geometry(data, profile$phi, matrix(1, 1L, 1L), factor = TRUE)
    ratios <- 2^(-(0:ceiling(2 * log2(data$n))) / 4)
    original <- sb_lasso_supports(factor, data$n, 1L, ratios)
    computed <- milk_curvature_supports(factor, data$n, 1L, ratios)
    if (!identical(original, computed)) stop("The actual-data curvature support check failed.")
    checks[[length(checks) + 1L]] <- data.table(tau = tau, h = h, columns = ncol(data$X),
      identical_supports = TRUE, paths = length(ratios))
  }
}
fwrite(rbindlist(checks), "results/preflight/milk_scs_20261031/curvature_solver_checks.csv")
print(rbindlist(checks))
