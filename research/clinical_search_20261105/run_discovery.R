library(data.table)
library(jsonlite)
library(yaml)
library(matrixStats)
library(glmnet)
library(quantreg)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 5L)
implementation <- normalizePath(args[5L], mustWork = TRUE)
config <- read_yaml(file.path(implementation, "config.yml"))
for (file in c("prediction_core.R", "ridge_direction.R", "panel.R")) source(file.path(implementation, file))
source("research/score_validity_20261026/load_joint.R")
source("research/joint_profile_manuscript_20261026/conic_search.R")
method <- args[2L]
if (startsWith(method, "POOLED")) source("research/soynam_application_20261031/load_pooled_lp.R")
outer <- file.path(config$output_root, args[1L])
rows <- fread(file.path(outer, "training.csv"))
evaluation <- fread(file.path(outer, "evaluation.csv"))
stopifnot(all(rows$role == "training"), all(evaluation$role == "discovery_evaluation"),
  !length(intersect(rows$strain_id, evaluation$strain_id)))
panel <- cp_panel(config$input_directory)
panel$recipes <- readRDS(file.path(outer, "recipes.rds"))
directory <- file.path(outer, "models", paste0(method, "_", args[3L]))
if (dir.exists(directory)) stop("拟合目录已经存在")
dir.create(directory, recursive = TRUE)
progress <- function(stage, details = list()) {
  record <- c(list(stage = stage, time = format(Sys.time(), tz = "UTC", usetz = TRUE)), details)
  write_json(record, file.path(directory, "progress.json"), auto_unbox = TRUE, pretty = TRUE)
  cat(stage, paste(unlist(details), collapse = " "), "\n")
  flush.console()
}
write_json(list(method = method, outer_fold = args[1L], config = config,
  code_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE), formal_validation_used = FALSE),
  file.path(directory, "request.json"), auto_unbox = TRUE, pretty = TRUE)
started <- proc.time()[["elapsed"]]
offset <- method %in% c("JOINT_OFFSET", "POOLED_OFFSET")
quantiles <- unlist(config$quantiles)
if (method == "PRESS_RIDGE") {
  direction <- panel$recipes$training$direction_recipe
  recipe <- direction$base_recipe
  features <- cp_features(evaluation, panel, recipe)
  oof <- direction$background$leave_out_prediction[rows$sample_id] * recipe$y_scale + recipe$y_center
  calibration <- lapply(quantiles, function(value) ep_calibrate(rows, oof, value, "strain_id", "condition"))
  base <- drop(ep_predict(features, direction$background$beta))
  prediction <- vapply(calibration, function(value) base + value[evaluation$condition], numeric(nrow(evaluation)))
  fit <- list(beta = direction$background$beta, recipe = recipe, calibration = calibration,
    prediction = prediction, quantiles = quantiles, background = direction$background)
} else if (method == "CLINICAL_QR") {
  recipe <- panel$recipes$training$direction_recipe$base_recipe
  features <- cp_features(rows, panel, recipe)
  test <- cp_features(evaluation, panel, recipe)
  clinical <- which(colnames(features$X) == "(Intercept)" |
    startsWith(colnames(features$X), "C_") | startsWith(colnames(features$X), "M_C_"))
  X <- features$X[, clinical, drop = FALSE]
  decomposition <- qr(X)
  keep <- sort(decomposition$pivot[seq_len(decomposition$rank)])
  X <- X[, keep, drop = FALSE]
  native <- lapply(quantiles, function(value) rq.fit(X * features$weight, features$y * features$weight, tau = value, method = "br"))
  beta <- do.call(cbind, lapply(native, `[[`, "coefficients"))
  prediction <- test$X[, colnames(X), drop = FALSE] %*% beta * recipe$y_scale + recipe$y_center
  fit <- list(beta = beta, recipe = recipe, prediction = prediction, quantiles = quantiles, native = native)
} else {
  recipe_function <- function(current, panel, config) cp_cached_recipe(current, panel, config, offset)
  solver <- if (method == "ENET") "ENET" else if (startsWith(method, "POOLED")) "POOLED" else "JOINT"
  fit <- ep_fit_prediction(rows, evaluation, panel, solver, as.numeric(args[4L]), config, progress, directory,
    recipe_function, if (offset) cp_offset_features else cp_features, "strain_id", "condition",
    if (offset) ep_predict_offset else ep_predict, if (method == "ENET") "oof" else "training")
}
stopifnot(all(is.finite(fit$prediction)))
fit$offset_model <- offset
saveRDS(fit, file.path(directory, "fit.rds"))
records <- rbindlist(lapply(seq_along(fit$quantiles), function(k)
  data.table(sample_id = evaluation$sample_id, strain_id = evaluation$strain_id,
    method = method, tau = fit$quantiles[k], prediction = fit$prediction[, k])))
fwrite(records, file.path(directory, "prediction.csv"))
writeLines(trimws(capture.output(warnings()), which = "right"), file.path(directory, "warnings.txt"))
write_json(list(status = "completed", stage = "discovery", runtime_sec = proc.time()[["elapsed"]] - started,
  formal_validation_used = FALSE, manuscript_adoption = FALSE), file.path(directory, "status.json"),
  auto_unbox = TRUE, pretty = TRUE)
