library(data.table)
library(jsonlite)
library(quantreg)
source("research/score_validity_20261026/load_joint.R")
source("research/soynam_application_20261031/load_pooled_lp.R")
source("research/drug_response_application_20261101/pooled.R")

objects <- new.env()
load("results/scratch/empirical_selection_20261101/fbioselect_source/test_gdsc_data.RData", objects)
fd <- objects$fd
basis <- splines::bs(fd$DosageSQ, df = 5L, degree = 3L, intercept = TRUE)
colnames(basis) <- paste0("(Intercept)_dose", 1:5)
X <- cbind(basis, as.matrix(fd[, 6:25]))
factor <- c(rep(0, 5L), rep(1, 20L))
data <- list(X = X, y = fd$Response, penalty_factor = factor)
lambda <- c(0.05, 0.025)
checks <- rbindlist(lapply(c(0.25, 0.5, 0.75), function(tau) {
  fit <- drug_pooled_path(data, tau, lambda = lambda)
  standardised <- sweep(sweep(X, 2L, fit$means, "-"), 2L, fit$scales, "/")
  rbindlist(lapply(seq_along(lambda), function(j) {
    reference <- rq.fit.lasso(standardised, fd$Response, tau = tau,
                               lambda = 2 * nrow(X) * lambda[j] * factor, eps = 1e-9)
    objective <- function(prediction, beta) {
      residual <- fd$Response - prediction
      mean(residual * (tau - (residual < 0))) + lambda[j] * sum(factor * abs(beta))
    }
    actual <- objective(drop(X %*% fit$beta[, j]), fit$beta[, j] * fit$scales)
    expected <- objective(drop(standardised %*% reference$coefficients), reference$coefficients)
    data.table(tau = tau, lambda = lambda[j], implementation_objective = actual,
      reference_objective = expected, objective_error = abs(actual - expected),
      gap = fit$diagnostics$gap[j], mapping_error = fit$diagnostics$mapping_error[j])
  }))
}))
stopifnot(max(checks$objective_error) < 1e-7, max(abs(checks$gap)) < 1e-7)
directory <- "results/preflight/empirical_selection_20261101/quantile_lasso"
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
fwrite(checks, file.path(directory, "objective_checks.csv"))
write_json(list(status = "passed", actual_observations_used = TRUE,
                maximum_objective_error = max(checks$objective_error),
                reference = "quantreg::rq.fit.lasso", implementation = "HiGHS_certified_pooled_path"),
           file.path(directory, "checks.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
print(checks)
