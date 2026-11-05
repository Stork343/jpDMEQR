library(yaml)
library(jsonlite)
source("research/clinical_quantile_20261105/engine.R")
source("research/clinical_quantile_20261105/orthogonal.R")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% c(0L, 2L, 3L))
config <- read_yaml(if (length(args) == 3L) args[1L] else "config/clinical_quantile_20261105.yml")
if (length(args) == 3L) args <- args[-1L]
study <- read_yaml(if (length(args)) args[1L] else "config/application_gnhs_baseline_bmi_20261105.yml")
if (length(args)) config$inference_target <- args[2L]
if (config$prediction_method == "CQ_SMOOTH") {
  reticulate::use_python(file.path(getwd(), "results/scratch/cq_solver_env/bin/python"), required = TRUE)
  Sys.setenv(NUMBA_CACHE_DIR = normalizePath("results/scratch/runtime", mustWork = TRUE))
  reticulate::source_python(file.path(config$output_root, "implementation", "smooth_path.py"))
}
base <- file.path(config$output_root, study$study_id)
root <- file.path(config$output_root, paste0("orthogonal_", study$study_id))
if (dir.exists(root)) stop("推断目录已经存在")
dir.create(root, recursive = TRUE)
panel <- cp_panel(study$input_directory)
tau <- config$inference_quantile
contributions <- diagnostics <- list()
for (fold in seq_len(study$outer_folds)) {
  outer <- file.path(study$output_root, paste0("outer_", fold))
  rows <- fread(file.path(outer, "training.csv"))
  evaluation <- fread(file.path(outer, "evaluation.csv"))
  stopifnot(!length(intersect(rows$strain_id, evaluation$strain_id)),
    all(evaluation$role == "discovery_evaluation"))
  central <- readRDS(file.path(base, paste0("outer_", fold), "tau075", "fit.rds"))
  training <- cq_features(rows, panel, central$recipe, central$columns)
  testing <- cq_features(evaluation, panel, central$recipe, central$columns)
  directory <- file.path(root, paste0("outer_", fold))
  dir.create(directory)
  delta <- min(config$density_spacing, 0.25 * training$n^(-1/5))
  fitted <- lapply(c(tau - delta, tau + delta), function(value) {
    lambda <- cq_lambda(training, value, central$alpha, config)
    index <- central$grid[central$choice]$index
    if (config$prediction_method == "CQ_ENET") return(cq_exact(training, value, lambda[index], central$alpha, config))
    stopifnot(config$prediction_method == "CQ_SMOOTH")
    fitted <- fit_smooth_path(training$X, training$y, training$weight, training$penalty_factor,
      c(lambda[index]), central$alpha, value, training$n, central$beta, config$solver_tolerance)
    beta <- setNames(fitted$beta[, 1L], colnames(training$X))
    list(beta = beta, tau = value, alpha = central$alpha, lambda = lambda[index],
      bandwidth = fitted$bandwidth, diagnostics = fitted$diagnostics)
  })
  density_train <- cq_density(training, fitted[[1L]]$beta, fitted[[2L]]$beta, delta, config)
  density_test <- cq_density(testing, fitted[[1L]]$beta, fitted[[2L]]$beta, delta, config)
  correlation <- cq_score_correlation(rows, central$oof, tau, config)
  saveRDS(list(density_fits = fitted, delta = delta, correlation = correlation), file.path(directory, "nuisance.rds"))
  for (label in c("independence", "exchangeable")) {
    rho <- if (label == "independence") 0 else correlation$working
    projection <- cq_projection(training, density_train$value, config$inference_target,
      rho, tau, rows$inner_fold, config)
    result <- cq_influence(testing, central$beta, density_test$value, projection)
    contributions[[length(contributions) + 1L]] <- data.table(fold, method = label,
      strain_id = result$patient, initial = result$initial, influence = result$influence)
    diagnostics[[length(diagnostics) + 1L]] <- data.table(fold, method = label, rho,
      raw_rho = correlation$raw, target = config$inference_target, tau, delta,
      sensitivity = projection$sensitivity, maximum_balance = projection$maximum_balance,
      clinical_balance = projection$clinical_balance, projection_lambda = projection$lambda,
      training_density_constrained = density_train$lower_bound_count + density_train$upper_bound_count,
      evaluation_density_constrained = density_test$lower_bound_count + density_test$upper_bound_count,
      training_rows = nrow(rows), evaluation_rows = nrow(evaluation))
    saveRDS(list(projection = projection, result = result), file.path(directory, paste0(label, ".rds")))
  }
  cat("orthogonal completed fold", fold, "\n")
  flush.console()
}
all <- rbindlist(contributions)
stopifnot(!anyDuplicated(all[, .(method, strain_id)]))
all[, centered := influence - mean(influence), by = .(method, fold)]
summary <- all[, .(patients = .N, estimate = mean(initial + influence),
  standard_error = sqrt(sum(centered^2)) / .N), by = method]
summary[, `:=`(lower = estimate - qnorm(0.975) * standard_error,
  upper = estimate + qnorm(0.975) * standard_error)]
fwrite(all, file.path(root, "patient_influence.csv"))
report <- config$report_directory
fwrite(summary, file.path(report, paste0(study$study_id, "_orthogonal_summary.csv")))
fwrite(rbindlist(diagnostics), file.path(report, paste0(study$study_id, "_orthogonal_diagnostics.csv")))
write_json(list(status = "completed", role = "exploratory_conditional_inference",
  target = config$inference_target, tau = tau, formal_validation_used = FALSE, manuscript_adoption = FALSE),
  file.path(root, "status.json"), auto_unbox = TRUE, pretty = TRUE)
print(summary)
