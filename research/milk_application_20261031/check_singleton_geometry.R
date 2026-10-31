library(data.table)
source("research/score_validity_20261026/load_joint.R")
cohort <- readRDS("results/raw/application/milk_scs_20261031/cohorts.rds")$eligible
counts <- table(cohort$cluster_id)
index <- which(cohort$cluster_id %in% names(counts)[counts == 1L])
data <- sb_data(cohort$y[index], cohort$X[index, 1:13], cohort$Z[index, , drop = FALSE],
  cohort$cluster_id[index])
if (any(lengths(data$rows) != 1L)) stop("The observed singleton selection is inconsistent.")
beta <- numeric(ncol(data$X))
beta[2:13] <- median(data$y)
direction <- seq_len(ncol(data$X)) / sqrt(sum(seq_len(ncol(data$X))^2))
checks <- list()
for (tau in c(0.25, 0.5, 0.75)) {
  for (lambda in c(0.5, 1, 2)) {
    for (h in length(unique(cohort$cluster_id))^(-0.3) * c(1, 0.5)) {
      Lambda <- matrix(lambda, 1L, 1L)
      profile <- ds_profile(data, beta, tau, h, Lambda)
      kappa <- profile$phi / (1 + profile$phi / lambda)
      explicit <- crossprod(data$X, data$X * as.numeric(kappa)) / data$n
      computed <- sb_geometry(data, profile$phi, Lambda)
      gradient <- function(coefficient) {
        current <- ds_profile(data, coefficient, tau, h, Lambda)
        -drop(crossprod(data$X, current$psi)) / data$n
      }
      derivative <- (gradient(beta + 1e-6 * direction) - gradient(beta - 1e-6 * direction)) / 2e-6
      errors <- c(hessian_identity = max(abs(explicit - computed)),
        finite_difference = max(abs(derivative - computed %*% direction)))
      tolerances <- c(1e-10, 2e-5)
      if (any(!is.finite(errors)) || any(errors > tolerances)) stop("The observed singleton geometry failed.")
      checks[[length(checks) + 1L]] <- data.table(singleton_cows = data$n, tau = tau,
        lambda = lambda, h = h, check = names(errors), error = as.numeric(errors), tolerance = tolerances)
    }
  }
}
fwrite(rbindlist(checks), "results/preflight/milk_scs_20261031/singleton_geometry_checks.csv")
print(rbindlist(checks)[, .(maximum_error = max(error)), by = check])
