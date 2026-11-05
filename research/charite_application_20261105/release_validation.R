library(data.table)
library(jsonlite)
library(yaml)
library(matrixStats)
source("research/charite_application_20261105/features.R")
source("research/charite_application_20261105/predict.R")
source("research/empirical_selection_20261101/prediction_core.R")
source("research/empirical_selection_20261101/ridge_direction.R")
config <- read_yaml("config/application_charite_hb_20261105.yml")
root <- config$output_root
release_path <- file.path(root, "validation_release.json")
if (file.exists(release_path)) stop("验证登记已经存在")
linear <- fread(file.path(root, "progress.csv"))
forest <- fread(file.path(root, "forest_queue", "progress.csv"))
stopifnot(nrow(linear) == 15L, nrow(forest) == 2L, all(linear$status == "completed"), all(forest$status == "completed"))
directories <- c(linear$directory, forest$directory, file.path(root, "models", c("XGB_RAW_all", "XGB_OFFSET_all")))
for (directory in directories) {
  status <- read_json(file.path(directory, "status.json"), simplifyVector = TRUE)
  stopifnot(status$status == "completed", status$stage == "training", !status$formal_validation_used)
}
for (filename in c("training_checks.json", "forest_checks.json", "boost_raw_checks.json", "boost_offset_checks.json")) {
  check <- read_json(file.path(config$report_directory, filename), simplifyVector = TRUE)
  stopifnot(check$status == "passed", !check$formal_validation_used)
}
rows <- fread(file.path(config$input_directory, "training.csv"))
panel <- ch_panel(config$input_directory)
slopes <- checks <- list()
for (value in unlist(config$quantiles)) {
  directory <- file.path(root, "models", sprintf("JOINT_OFFSET_tau%03d", as.integer(value * 100)))
  fit <- readRDS(file.path(directory, "fit.rds"))
  expanded <- ch_expand_offset(fit$beta, fit$recipe, fit$calibration[[1L]]["hemoglobin"])
  base <- fit$recipe$direction_recipe$base_recipe
  genes <- base$genes[paste0("P_", base$genes) %in% names(expanded)]
  columns <- paste0("P_", genes)
  slopes[[length(slopes) + 1L]] <- data.table(protein = genes, tau = value,
    total_slope = expanded[columns] / base$scale[columns], observed_fraction = base$observed[genes])
  inner <- readRDS(file.path(directory, "inner_1.rds"))
  fitting <- rows[inner_fold != 1L]
  held <- rows[inner_fold == 1L]
  fitting_features <- ch_offset_features(fitting, panel, inner$recipe)
  beta <- setNames(as.numeric(inner$path$beta[, fit$choice]), colnames(fitting_features$X))
  calibration <- ep_calibrate(fitting, drop(ep_predict_offset(fitting_features, beta)), value, "strain_id", "condition")
  held_features <- ch_features(held, panel, inner$recipe$direction_recipe$base_recipe)
  inner_expanded <- ch_expand_offset(beta, inner$recipe, calibration["hemoglobin"])
  algebraic <- drop(held_features$X %*% inner_expanded)
  predictor <- ch_predict(list(beta = beta, recipe = inner$recipe, calibration = list(calibration),
    quantiles = value, offset_model = TRUE, model_type = "linear"), held[, !"y"], panel)
  error <- max(abs(algebraic - inner$prediction[, fit$choice]), abs(predictor[, 1L] - inner$prediction[, fit$choice]))
  stopifnot(error < 1e-9)
  checks[[length(checks) + 1L]] <- data.table(tau = value, internal_heldout_patients = uniqueN(held$strain_id), maximum_error = error)
}
slopes <- rbindlist(slopes)
wide <- dcast(slopes, protein + observed_fraction ~ tau, value.var = "total_slope")
wide[, difference := `0.25` - `0.75`]
prespecified <- unlist(config$inference$prespecified_proteins)
eligible <- wide[!protein %in% prespecified & !grepl(";", protein, fixed = TRUE) &
  observed_fraction >= config$inference$selection_minimum_observed_fraction]
eligible[, magnitude := abs(difference)]
setorder(eligible, -magnitude, protein)
additional <- head(eligible$protein, config$inference$additional_training_selected_proteins)
stopifnot(length(additional) == config$inference$additional_training_selected_proteins,
  all(prespecified %in% wide$protein))
targets <- data.table(protein = c(prespecified, additional),
  selection = c(rep("prespecified", length(prespecified)), rep("training_selected", length(additional))))
targets <- merge(targets, wide[, .(protein, training_slope_difference = difference)], by = "protein", sort = FALSE)
fwrite(targets, file.path(config$report_directory, "frozen_inference_targets.csv"))
fwrite(slopes, file.path(config$report_directory, "training_total_slopes.csv"))
fwrite(rbindlist(checks), file.path(config$report_directory, "coefficient_expansion_checks.csv"))
ledger_path <- "results/preflight/empirical_selection_20261101/validation_attempts.csv"
if (file.exists(ledger_path)) stop("验证清单需要核查")
fwrite(data.table(attempt_index = config$validation$attempt_index, study_id = config$study_id,
  outcome = config$response, decision_quantile = config$decision_quantile, total_alpha = config$validation$total_alpha,
  prediction_alpha = config$validation$prediction_alpha, inference_alpha = config$validation$inference_alpha,
  status = "registered", registered_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)), ledger_path)
write_json(list(status = "frozen", study_id = config$study_id, training_patients = uniqueN(rows$strain_id),
  model_directories = directories, target_proteins = targets$protein, primary_quantile = config$decision_quantile,
  validation = config$validation, primary_prediction_rule = config$primary_prediction_rule,
  code_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  frozen_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), formal_validation_used = FALSE),
  release_path, auto_unbox = TRUE, pretty = TRUE, digits = NA)
print(targets)
