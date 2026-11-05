library(data.table)
library(jsonlite)
library(yaml)
library(grf)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 5L)
fold_name <- args[1L]
method <- args[2L]
implementation <- normalizePath(args[5L], mustWork = TRUE)
config <- read_yaml(file.path(implementation, "config.yml"))
for (file in c("prediction_core.R", "ridge_direction.R", "features.R", "forest.R")) source(file.path(implementation, file))
stopifnot(method %in% config$methods)
fold <- as.integer(sub("fold", "", fold_name))
panel <- ibd_panel(config$input_directory)
rows <- fread(file.path(config$input_directory, "discovery.csv"))
training <- rows[outer_fold != fold]
evaluation <- rows[outer_fold == fold]
roles <- fread(file.path(config$source_root, fold_name, "training_roles.csv"))
stopifnot(identical(training$sample_id, roles$sample_id), all(rows$role == "discovery"))
training[, inner_fold := roles$inner_fold]
recipes <- readRDS(file.path(config$source_root, fold_name, "recipes.rds"))
directory <- file.path(config$output_root, fold_name, paste0(method, "_native"))
if (dir.exists(directory)) stop("运行目录已经存在")
dir.create(directory, recursive = TRUE)
progress <- function(stage, details = list()) {
  record <- c(list(stage = stage, time = format(Sys.time(), tz = "UTC", usetz = TRUE)), details)
  write_json(record, file.path(directory, "progress.json"), auto_unbox = TRUE, pretty = TRUE)
  cat(record$time, stage, paste(unlist(details), collapse = " "), "\n")
  flush.console()
}
write_json(list(method = method, fold = fold, config = config, phase = "method_development",
  implementation_directory = implementation, code_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  formal_validation_used = FALSE, manuscript_adoption = FALSE), file.path(directory, "request.json"),
  auto_unbox = TRUE, pretty = TRUE)
started <- proc.time()[["elapsed"]]
fit <- rs_forest_prediction(training, evaluation, panel, recipes, method == "QRF_OFFSET", progress,
  directory, base_feature_function = ibd_features)
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
