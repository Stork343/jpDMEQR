library(data.table)
library(jsonlite)
library(yaml)
library(glmnet)

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 5L)
assay_name <- args[1L]
method <- args[2L]
run_id <- args[3L]
implementation <- normalizePath(args[5L], mustWork = TRUE)
config <- read_yaml(file.path(implementation, "config.yml"))
base <- read_yaml(file.path(implementation, "base_config.yml"))
stopifnot(assay_name %in% names(base$assays))
base$conditions <- base$assays[[assay_name]]
stopifnot(method %in% c("JOINT_DIRECTION", "POOLED_DIRECTION", "PRESS_RIDGE", "JOINT_OFFSET", "POOLED_OFFSET"))
tau <- if (method == "PRESS_RIDGE") NA_real_ else as.numeric(args[4L])
source(file.path(implementation, "prediction_core.R"))
source(file.path(implementation, "ridge_direction.R"))
source(file.path(implementation, "contrasts.R"))
source(file.path(implementation, "features.R"))
source("research/score_validity_20261026/load_joint.R")
source("research/joint_profile_manuscript_20261026/conic_search.R")
if (startsWith(method, "POOLED_")) source("research/soynam_application_20261031/load_pooled_lp.R")
panel <- proteome_panel(config$input_directory)
recipe_path <- if (isTRUE(config$recipes_by_assay)) file.path(config$output_root, assay_name, "recipes.rds") else
  file.path(config$output_root, "recipes.rds")
panel$direction_recipes <- readRDS(recipe_path)
rows <- fread(file.path(config$input_directory, assay_name, "discovery.csv"))
training <- rows[exploration_role == "training"]
evaluation <- rows[exploration_role == "evaluation"]
stopifnot(all(rows$role == "discovery"), !length(intersect(training$strain_id, evaluation$strain_id)))
directory <- file.path(config$output_root, assay_name, paste0(method, "_", run_id))
if (dir.exists(directory)) stop("运行目录已经存在")
dir.create(directory, recursive = TRUE)
progress <- function(stage, details = list()) {
  record <- c(list(stage = stage, method = method, time = format(Sys.time(), tz = "UTC", usetz = TRUE)), details)
  write_json(record, file.path(directory, "progress.json"), auto_unbox = TRUE, pretty = TRUE)
  cat(record$time, stage, paste(unlist(details), collapse = " "), "\n")
  flush.console()
}
write_json(list(method = method, tau = tau, phase = "method_development", config = config,
                implementation_directory = implementation, code_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
                formal_validation_used = FALSE, manuscript_adoption = FALSE),
           file.path(directory, "request.json"), auto_unbox = TRUE, pretty = TRUE, na = "null")
recipe_function <- function(rows, panel, config) {
  ids <- sort(unique(rows$strain_id))
  keys <- names(panel$direction_recipes)[vapply(panel$direction_recipes,
    function(recipe) identical(recipe$training_strains, ids), logical(1))]
  stopifnot(length(keys) == 1L)
  panel$direction_recipes[[keys]]
}
feature_function <- if (endsWith(method, "OFFSET")) {
  function(rows, panel, recipe) ep_offset_features(rows, panel, recipe, proteome_features)
} else {
  function(rows, panel, recipe) ep_direction_features(rows, panel, recipe, proteome_features)
}
started <- proc.time()[["elapsed"]]
if (method == "PRESS_RIDGE") {
  recipe <- panel$direction_recipes$training
  if (!is.null(recipe$direction_recipe)) recipe <- recipe$direction_recipe
  background <- recipe$background
  quantiles <- unlist(config$quantiles)
  oof <- background$leave_out_prediction[training$sample_id] * recipe$y_scale + recipe$y_center
  calibration <- lapply(quantiles, function(value) ep_calibrate(training, oof, value, "strain_id", "condition"))
  test <- proteome_features(evaluation, panel, recipe$base_recipe)
  predicted <- drop(ep_predict(test, background$beta))
  prediction <- vapply(calibration, function(offset) predicted + offset[evaluation$condition], numeric(nrow(evaluation)))
  fit <- list(beta = background$beta, recipe = recipe$base_recipe, prediction = prediction, quantiles = quantiles,
              calibration = calibration, oof = oof, choice = background$choice,
              grid = data.table(lambda = background$lambda), loss = background$loss)
} else {
  solver_method <- if (startsWith(method, "JOINT_")) "JOINT" else "POOLED"
  prediction_function <- if (endsWith(method, "OFFSET")) ep_predict_offset else ep_predict
  calibration_mode <- if (isTRUE(config$calibrated_selection)) "training" else "oof"
  fit <- ep_fit_prediction(training, evaluation, panel, solver_method, tau, base, progress, directory,
                           recipe_function, feature_function, "strain_id", "condition", prediction_function, calibration_mode)
}
elapsed <- proc.time()[["elapsed"]] - started
saveRDS(fit, file.path(directory, "fit.rds"))
predictions <- rbindlist(lapply(seq_along(fit$quantiles), function(j)
  data.table(sample_id = evaluation$sample_id, strain_id = evaluation$strain_id, condition = evaluation$condition,
    y = evaluation$y, tau = fit$quantiles[j], prediction = fit$prediction[, j],
    loss = ep_loss(evaluation$y, fit$prediction[, j], fit$quantiles[j]))))
fwrite(predictions, file.path(directory, "predictions.csv"))
losses <- predictions[, .(loss = mean(loss)), by = .(strain_id, tau)]
fwrite(losses, file.path(directory, "strain_loss.csv"))
summary <- losses[, .(strains = .N, loss = mean(loss), se = sd(loss) / sqrt(.N)), by = tau]
summary[, `:=`(method = method, runtime_sec = elapsed, phase = "method_development")]
fwrite(summary, file.path(directory, "summary.csv"))
write_json(list(status = "completed", runtime_sec = elapsed, formal_validation_used = FALSE, manuscript_adoption = FALSE),
           file.path(directory, "status.json"), auto_unbox = TRUE, pretty = TRUE)
progress("completed", list(runtime_sec = elapsed))
print(summary)
