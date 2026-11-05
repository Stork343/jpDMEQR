library(data.table)
library(jsonlite)
library(yaml)
library(glmnet)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 5L)
fold_name <- args[1L]
method <- args[2L]
run_id <- args[3L]
tau <- as.numeric(args[4L])
implementation <- normalizePath(args[5L], mustWork = TRUE)
config <- read_yaml(file.path(implementation, "config.yml"))
for (file in c("prediction_core.R", "ridge_direction.R", "features.R")) source(file.path(implementation, file))
source("research/score_validity_20261026/load_joint.R")
source("research/joint_profile_manuscript_20261026/conic_search.R")
if (method == "POOLED_OFFSET") source("research/soynam_application_20261031/load_pooled_lp.R")
stopifnot(method %in% c("JOINT_OFFSET", "PRESS_RIDGE", "ENET", "POOLED_OFFSET", "JOINT_DIRECT"))
fold <- as.integer(sub("fold", "", fold_name))
panel <- ibd_panel(config$input_directory)
panel$recipes <- readRDS(file.path(config$output_root, fold_name, "recipes.rds"))
rows <- fread(file.path(config$input_directory, "discovery.csv"))
training <- rows[outer_fold != fold]
evaluation <- rows[outer_fold == fold]
roles <- fread(file.path(config$output_root, fold_name, "training_roles.csv"))
stopifnot(identical(roles$sample_id, training$sample_id), all(rows$role == "discovery"))
training[, inner_fold := roles$inner_fold]
directory <- file.path(config$output_root, fold_name, paste0(method, "_", run_id))
if (dir.exists(directory)) stop("运行目录已经存在")
dir.create(directory, recursive = TRUE)
progress <- function(stage, details = list()) {
  record <- c(list(stage = stage, time = format(Sys.time(), tz = "UTC", usetz = TRUE)), details)
  write_json(record, file.path(directory, "progress.json"), auto_unbox = TRUE, pretty = TRUE)
  cat(record$time, stage, paste(unlist(details), collapse = " "), "\n")
  flush.console()
}
write_json(list(method = method, tau = tau, fold = fold, config = config,
  code_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE), implementation_directory = implementation,
  formal_validation_used = FALSE, manuscript_adoption = FALSE), file.path(directory, "request.json"),
  auto_unbox = TRUE, pretty = TRUE, na = "null")
recipe_function <- function(current, panel, config) {
  keys <- names(panel$recipes)[vapply(panel$recipes,
    function(recipe) identical(recipe$training_strains, sort(unique(current$strain_id))), logical(1))]
  stopifnot(length(keys) == 1L)
  recipe <- panel$recipes[[keys]]
  if (method %in% c("ENET", "JOINT_DIRECT")) recipe$direction_recipe$base_recipe else recipe
}
offset_features <- function(current, panel, recipe) ep_offset_features(current, panel, recipe, ibd_features)
started <- proc.time()[["elapsed"]]
if (method == "PRESS_RIDGE") {
  direction <- panel$recipes$training$direction_recipe
  background <- direction$background
  quantiles <- unlist(config$quantiles)
  oof <- background$leave_out_prediction[training$sample_id] * direction$y_scale + direction$y_center
  calibration <- lapply(quantiles, function(value) ep_calibrate(training, oof, value, "strain_id", "condition"))
  test <- ibd_features(evaluation, panel, direction$base_recipe)
  base <- drop(ep_predict(test, background$beta))
  prediction <- vapply(calibration, function(value) base + value[evaluation$condition], numeric(nrow(evaluation)))
  fit <- list(beta = background$beta, recipe = direction$base_recipe, prediction = prediction,
    quantiles = quantiles, calibration = calibration, oof = oof, choice = background$choice,
    grid = data.table(lambda = background$lambda), loss = background$loss)
} else {
  offset <- method %in% c("JOINT_OFFSET", "POOLED_OFFSET")
  fit <- ep_fit_prediction(training, evaluation, panel,
    if (method == "ENET") "ENET" else if (method == "POOLED_OFFSET") "POOLED" else "JOINT",
    tau, config, progress, directory, recipe_function, if (offset) offset_features else ibd_features,
    "strain_id", "condition", if (offset) ep_predict_offset else ep_predict,
    if (method == "ENET") "oof" else "training")
}
elapsed <- proc.time()[["elapsed"]] - started
saveRDS(fit, file.path(directory, "fit.rds"))
predictions <- rbindlist(lapply(seq_along(fit$quantiles), function(k)
  data.table(sample_id = evaluation$sample_id, strain_id = evaluation$strain_id, condition = evaluation$condition,
    y = evaluation$y, tau = fit$quantiles[k], prediction = fit$prediction[, k],
    loss = ep_loss(evaluation$y, fit$prediction[, k], fit$quantiles[k]))))
fwrite(predictions, file.path(directory, "predictions.csv"))
summary <- predictions[, .(loss = mean(loss)), by = .(strain_id, tau)][,
  .(participants = .N, loss = mean(loss), se = sd(loss) / sqrt(.N)), by = tau]
summary[, `:=`(method = method, fold = fold, runtime_sec = elapsed)]
fwrite(summary, file.path(directory, "summary.csv"))
capture.output(sessionInfo(), file = file.path(directory, "session.txt"))
write_json(list(status = "completed", runtime_sec = elapsed, formal_validation_used = FALSE,
  manuscript_adoption = FALSE), file.path(directory, "status.json"), auto_unbox = TRUE, pretty = TRUE)
progress("completed", list(runtime_sec = elapsed))
