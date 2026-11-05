library(yaml)
library(jsonlite)
source("research/clinical_quantile_20261105/engine.R")
args <- commandArgs(trailingOnly = TRUE)
config <- read_yaml(if (length(args)) args[1L] else "config/clinical_quantile_20261105.yml")
output <- config$report_directory
aggregate <- fread(file.path(output, "adaptive_inference.csv"))
checks <- list()
for (j in seq_len(nrow(aggregate))) {
  record <- aggregate[j]
  root <- file.path(config$output_root, paste0("orthogonal_", record$study_id))
  patients <- fread(file.path(root, "patient_influence.csv"))
  first <- patients[method == "independence"][order(strain_id)]
  second <- patients[method == "exchangeable"][order(strain_id)]
  stopifnot(identical(first$strain_id, second$strain_id), identical(first$fold, second$fold))
  weight <- record$exchangeable_weight
  estimate <- (1 - weight) * mean(first$initial + first$influence) +
    weight * mean(second$initial + second$influence)
  se <- sqrt(sum(((1 - weight) * first$centered + weight * second$centered)^2)) / nrow(first)
  endpoint <- c(sqrt(sum(first$centered^2)), sqrt(sum(second$centered^2))) / nrow(first)
  stopifnot(abs(estimate - record$estimate) < 1e-12, abs(se - record$standard_error) < 1e-12,
    se <= min(endpoint) + 1e-12)
  checks[[length(checks) + 1L]] <- data.table(study_id = record$study_id,
    estimate_error = abs(estimate - record$estimate), standard_error_difference = abs(se - record$standard_error))
}
study <- read_yaml("config/application_ra_westchina_20261105.yml")
panel <- cp_panel(study$input_directory)
lambda_errors <- numeric()
for (fold in seq_len(study$outer_folds)) {
  rows <- fread(file.path(study$output_root, paste0("outer_", fold), "training.csv"))
  central <- readRDS(file.path(config$output_root, study$study_id, paste0("outer_", fold),
    sprintf("tau%03d", as.integer(100 * config$inference_quantile)), "fit.rds"))
  features <- cq_features(rows, panel, central$recipe, central$columns)
  nuisance <- readRDS(file.path(config$output_root, paste0("orthogonal_", study$study_id),
    paste0("outer_", fold), "nuisance.rds"))
  for (j in seq_along(nuisance$density_fits)) {
    fit <- nuisance$density_fits[[j]]
    lambda <- cq_lambda(features, fit$tau, fit$alpha, config)
    error <- abs(lambda[central$grid[central$choice]$index] - fit$lambda)
    stopifnot(error < 1e-12)
    lambda_errors <- c(lambda_errors, error)
  }
}
fwrite(rbindlist(checks), file.path(output, "independent_inference_checks.csv"))
write_json(list(status = "passed", density_direct_fit_comparisons = length(lambda_errors),
  maximum_density_lambda_difference = max(lambda_errors),
  maximum_aggregate_estimate_error = max(rbindlist(checks)$estimate_error),
  maximum_aggregate_se_error = max(rbindlist(checks)$standard_error_difference)),
  file.path(output, "inference_verification.json"), auto_unbox = TRUE, pretty = TRUE)
print(rbindlist(checks))
