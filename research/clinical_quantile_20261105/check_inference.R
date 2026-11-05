library(yaml)
library(jsonlite)
source("research/clinical_quantile_20261105/engine.R")
source("research/clinical_quantile_20261105/orthogonal.R")
config <- read_yaml("config/clinical_quantile_20261105.yml")
study <- read_yaml("config/application_gnhs_baseline_bmi_20261105.yml")
outer <- file.path(study$output_root, "outer_1")
rows <- fread(file.path(outer, "training.csv"))
evaluation <- fread(file.path(outer, "evaluation.csv"))
panel <- cp_panel(study$input_directory)
recipe <- readRDS(file.path(outer, "recipes.rds"))$training$direction_recipe$base_recipe
training <- cq_features(rows, panel, recipe)
columns <- colnames(training$X)[c(which(training$penalty_factor == 0),
  which(training$penalty_factor == 1)[seq_len(12L)])]
training <- cq_features(rows, panel, recipe, columns)
testing <- cq_features(evaluation, panel, recipe, columns)
tau <- config$inference_quantile
delta <- min(config$density_spacing, 0.25 * training$n^(-1/5))
fitted <- lapply(c(tau - delta, tau, tau + delta), function(value)
  cq_exact(training, value, 0.01, 0.5, config))
density_train <- cq_density(training, fitted[[1L]]$beta, fitted[[3L]]$beta, delta, config)
density_test <- cq_density(testing, fitted[[1L]]$beta, fitted[[3L]]$beta, delta, config)
records <- list()
for (rho in c(0, 0.4)) {
  projection <- cq_projection(training, density_train$value, config$inference_target,
    rho, tau, rows$inner_fold, config)
  result <- cq_influence(testing, fitted[[2L]]$beta, density_test$value, projection)
  stopifnot(length(result$influence) == uniqueN(evaluation$strain_id),
    all(is.finite(result$influence)), is.finite(result$initial),
    setequal(result$patient, evaluation$strain_id))
  records[[length(records) + 1L]] <- data.table(rho,
    training_patients = training$n, evaluation_patients = testing$n,
    projection_lambda = projection$lambda, sensitivity = projection$sensitivity,
    maximum_balance = projection$maximum_balance, clinical_balance = projection$clinical_balance,
    influence_sd = sd(result$influence), finite_contributions = length(result$influence))
}
output <- "results/preflight/clinical_quantile_20261105"
fwrite(rbindlist(records), file.path(output, "inference_implementation_checks.csv"))
write_json(list(status = "passed", role = "implementation_test_on_actual_discovery_patients",
  selected_test_proteins = 12L, target = config$inference_target, tau = tau,
  formal_validation_used = FALSE, manuscript_adoption = FALSE),
  file.path(output, "inference_implementation_acceptance.json"), auto_unbox = TRUE, pretty = TRUE)
print(rbindlist(records))
