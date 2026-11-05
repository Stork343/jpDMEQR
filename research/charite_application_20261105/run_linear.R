library(data.table)
library(jsonlite)
library(yaml)
library(matrixStats)
library(glmnet)
library(quantreg)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 5L, args[1L] == "hemoglobin")
method <- args[2L]
run_id <- args[3L]
tau <- as.numeric(args[4L])
implementation <- normalizePath(args[5L], mustWork = TRUE)
config <- read_yaml(file.path(implementation, "config.yml"))
for (file in c("prediction_core.R", "ridge_direction.R", "features.R", "train_only.R")) source(file.path(implementation, file))
source("research/score_validity_20261026/load_joint.R")
source("research/joint_profile_manuscript_20261026/conic_search.R")
if (method %in% c("POOLED", "POOLED_OFFSET")) source("research/soynam_application_20261031/load_pooled_lp.R")
rows <- fread(file.path(config$input_directory, "training.csv"))
stopifnot(all(rows$role == "training"))
panel <- ch_panel(config$input_directory)
panel$recipes <- readRDS(file.path(config$output_root, "recipes.rds"))
directory <- file.path(config$output_root, "models", paste0(method, "_", run_id))
if (dir.exists(directory)) stop("运行目录已经存在")
dir.create(directory, recursive = TRUE)
progress <- function(stage, details = list()) {
  record <- c(list(stage = stage, time = format(Sys.time(), tz = "UTC", usetz = TRUE)), details)
  write_json(record, file.path(directory, "progress.json"), auto_unbox = TRUE, pretty = TRUE)
  cat(record$time, stage, paste(unlist(details), collapse = " "), "\n")
  flush.console()
}
write_json(list(method = method, tau = tau, stage = "training", config = config,
  code_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  implementation_directory = implementation, formal_validation_used = FALSE),
  file.path(directory, "request.json"), auto_unbox = TRUE, pretty = TRUE, na = "null")
started <- proc.time()[["elapsed"]]
offset <- method %in% c("JOINT_OFFSET", "POOLED_OFFSET")
recipe_function <- function(current, panel, config) ch_cached_recipe(current, panel, config, offset)
if (method == "PRESS_RIDGE") {
  direction <- panel$recipes$training$direction_recipe
  background <- direction$background
  quantiles <- unlist(config$quantiles)
  oof <- background$leave_out_prediction[rows$sample_id] * direction$y_scale + direction$y_center
  calibration <- lapply(quantiles, function(value) ep_calibrate(rows, oof, value, "strain_id", "condition"))
  fit <- list(beta = background$beta, recipe = direction$base_recipe, quantiles = quantiles,
    calibration = calibration, choice = background$choice, grid = data.table(lambda = background$lambda),
    loss = background$loss, model_type = "linear")
} else if (method == "CLINICAL_QR") {
  recipe <- panel$recipes$training$direction_recipe$base_recipe
  features <- ch_features(rows, panel, recipe)
  clinical <- which(colnames(features$X) == "(Intercept)" | startsWith(colnames(features$X), "C_"))
  C <- features$X[, clinical, drop = FALSE]
  decomposition <- qr(C)
  keep <- sort(decomposition$pivot[seq_len(decomposition$rank)])
  C <- C[, keep, drop = FALSE]
  quantiles <- unlist(config$quantiles)
  native <- lapply(quantiles, function(value) rq.fit(C * features$weight, features$y * features$weight, tau = value, method = "br"))
  beta <- do.call(cbind, lapply(native, `[[`, "coefficients"))
  rownames(beta) <- colnames(C)
  stopifnot(all(is.finite(beta)))
  fit <- list(beta = beta, recipe = recipe, quantiles = quantiles, clinical_columns = colnames(C),
    native = native, model_type = "clinical", calibration = lapply(quantiles, function(value) c(hemoglobin = 0)))
} else {
  solver_method <- if (method == "ENET") "ENET" else if (startsWith(method, "POOLED")) "POOLED" else "JOINT"
  fit <- ch_train_prediction(rows, panel, solver_method, tau, config, progress, directory,
    recipe_function, if (offset) ch_offset_features else ch_features,
    if (offset) ep_predict_offset else ep_predict, if (method == "ENET") "oof" else "training")
  fit$model_type <- "linear"
}
fit$offset_model <- offset
elapsed <- proc.time()[["elapsed"]] - started
saveRDS(fit, file.path(directory, "fit.rds"))
capture.output(sessionInfo(), file = file.path(directory, "session.txt"))
capture.output(warnings(), file = file.path(directory, "warnings.txt"))
write_json(list(status = "completed", stage = "training", runtime_sec = elapsed,
  formal_validation_used = FALSE, manuscript_adoption = FALSE), file.path(directory, "status.json"),
  auto_unbox = TRUE, pretty = TRUE)
progress("training_completed", list(runtime_sec = elapsed))
