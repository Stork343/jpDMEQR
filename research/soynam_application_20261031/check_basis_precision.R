library(data.table)
library(yaml)
source("research/soynam_application_20261031/engine.R")
config <- read_yaml("config/application_soynam_20261031.yml")
root <- file.path(config$results_directory, "inference", config$execution_version)
checks <- list()
for (case_id in c("reference_q50", "bandwidth_075_q50")) {
  prepared <- readRDS(file.path(root, case_id, "prepared.rds"))
  request <- readRDS(file.path(root, case_id, "request.rds"))
  tau <- request$job$tau
  h <- request$job$h_multiplier * prepared$n^(-.3)
  reference <- if (case_id == "reference_q50") readRDS(file.path(root, case_id, "fit.rds")) else NULL
  for (fold in seq_along(prepared$folds)) for (coordinate in seq_along(config$coordinates)) {
    model <- prepared$folds[[fold]]$models[[coordinate]]
    fit <- ds_pair_small(model$data, model$initial, tau, h,
      matrix(request$job$lambda, 1L, 1L), gram = model$gram)
    for (bandwidth in c("full", "half")) {
      current <- fit[[bandwidth]]
      if (current$kkt > 1e-8 || current$prediction_error > 1e-7)
        stop("The refined coefficient system fails its original-coordinate checks.")
      coefficient_difference <- objective_difference <- influence_difference <- NA_real_
      if (!is.null(reference)) {
        previous <- reference$details[[coordinate]]$fits[[fold]][[bandwidth]]
        coefficient_difference <- max(abs(current$beta - previous$beta))
        objective_difference <- abs(current$objective - previous$objective)
        influence_difference <- max(abs(current$projected_scores - previous$projected_scores))
        if (coefficient_difference > 1e-5 || objective_difference > 1e-8 || influence_difference > 1e-5)
          stop("The precision calculation changed the fitted reference beyond numerical tolerance.")
      }
      checks[[length(checks) + 1L]] <- data.table(case_id, fold, coordinate, bandwidth,
        original_kkt = current$kkt, basis_kkt = current$basis_kkt,
        prediction_error = current$prediction_error, iterations = current$iterations,
        coefficient_difference, objective_difference, influence_difference)
    }
  }
}
fwrite(rbindlist(checks), "results/preflight/soynam_20261031/basis_precision_checks.csv")
print(rbindlist(checks))
