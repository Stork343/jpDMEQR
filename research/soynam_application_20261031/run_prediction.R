library(data.table)
library(yaml)
library(hqreg)
config <- read_yaml("config/application_soynam_20261031.yml")
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L || !arguments[1L] %in% unlist(config$prediction$methods))
  stop("Provide a registered prediction method and quantile.")
method <- arguments[1L]
tau <- as.numeric(arguments[2L])
if (!tau %in% unlist(config$quantiles)) stop("The requested prediction quantile is not registered.")
case_id <- paste0(method, "_q", sprintf("%02d", 100 * tau))
directory <- file.path(config$results_directory, "prediction", config$execution_version, case_id)
if (dir.exists(directory)) stop("The prediction case already exists.")
dir.create(directory, recursive = TRUE)
backup <- file.path(config$backup_directory, "prediction", config$execution_version, case_id)
dir.create(backup, recursive = TRUE, showWarnings = FALSE)
started <- Sys.time()
options(error = function() {
  saveRDS(list(case_id = case_id, error = geterrmessage(), calls = lapply(sys.calls(), deparse),
    started = as.character(started), failed = as.character(Sys.time())), file.path(directory, "failure.rds"))
  file.copy(file.path(directory, "failure.rds"), backup, overwrite = FALSE)
  quit(save = "no", status = 1L, runLast = FALSE)
})
source("research/soynam_application_20261031/engine.R")
source("research/soynam_application_20261031/prediction.R")
soy_context <- list(case_id = case_id, directory = directory)
input <- readRDS(file.path(config$results_directory, "input.rds"))
data <- soy_design(input, prediction = TRUE, canopy = method != "joint_genomic")
heldout <- input$line_folds[prediction_fold == config$prediction$heldout_fold, cluster_id]
training_rows <- which(!data$cluster_id %in% heldout)
test_rows <- which(data$cluster_id %in% heldout)
outer <- soy_subset(data, training_rows)
inner_assignment <- pc_folds(outer$rows, config$prediction$inner_folds)
inner <- lapply(seq_len(config$prediction$inner_folds), function(fold) {
  list(train = soy_subset(outer, unlist(outer$rows[inner_assignment != fold], use.names = FALSE)),
    validation = soy_subset(outer, unlist(outer$rows[inner_assignment == fold], use.names = FALSE)))
})
if (length(intersect(unique(outer$cluster_id), heldout))) stop("Prediction sample roles overlap.")
saveRDS(list(case_id = case_id, method = method, tau = tau, config = config,
  training_lines = unique(outer$cluster_id), test_lines = heldout,
  source_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE), started = as.character(started)),
  file.path(directory, "request.rds"))
fitted <- if (method == "pooled_quantile_lasso") {
  soy_pooled_prediction(outer, inner, tau, config$prediction)
} else {
  soy_joint_prediction(outer, inner, tau, config$prediction$working_lambda)
}
training_years <- data$observations$Year[training_rows]
calibration <- setNames(vapply(unlist(config$years), function(year) {
  stats::quantile((outer$y - fitted$oof)[training_years == year], tau, type = 1, names = FALSE)
}, numeric(1)), unlist(config$years))
fitted$calibration <- calibration
saveRDS(fitted, file.path(directory, "model.rds"), compress = FALSE)
raw_prediction <- drop(data$X[test_rows, , drop = FALSE] %*% fitted$beta)
test_metadata <- data$observations[test_rows]
prediction <- raw_prediction + calibration[as.character(test_metadata$Year)]
predictions <- data.table(cluster_id = data$cluster_id[test_rows], ID = test_metadata$ID,
  Year = test_metadata$Year, method = method, tau = tau,
  observed = data$y[test_rows] * config$response_divisor,
  raw_prediction = raw_prediction * config$response_divisor,
  prediction = prediction * config$response_divisor)
if (any(!is.finite(predictions$prediction))) stop("Nonfinite heldout predictions.")
predictions[, `:=`(loss = (observed - prediction) * (tau - (observed < prediction)),
  raw_loss = (observed - raw_prediction) * (tau - (observed < raw_prediction)),
  covered = as.integer(observed <= prediction), raw_covered = as.integer(observed <= raw_prediction))]
fwrite(predictions, file.path(directory, "predictions.csv"))
fwrite(predictions[, .(mean_loss = mean(loss), raw_mean_loss = mean(raw_loss),
  calibration = mean(covered), raw_calibration = mean(raw_covered)), by = .(Year, method, tau)],
  file.path(directory, "metrics.csv"))
jsonlite::write_json(list(case_id = case_id, status = "complete", method = method, tau = tau,
  training_lines = outer$n, test_lines = length(heldout),
  wall_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")), completed = as.character(Sys.time())),
  file.path(directory, "completion.json"), auto_unbox = TRUE, pretty = TRUE)
if (!all(file.copy(list.files(directory, full.names = TRUE), backup, overwrite = FALSE))) stop("Prediction backup failed.")
cat(case_id, "completed\n")
