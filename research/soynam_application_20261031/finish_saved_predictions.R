library(data.table)
library(jsonlite)
source("research/soynam_application_20261031/data.R")
root <- "results/raw/application/soynam_yield_20261031/prediction/paired_qp_20261031"
for (case_id in c("joint_canopy_q75", "joint_genomic_q75")) {
  directory <- file.path(root, case_id)
  if (file.exists(file.path(directory, "completion.json"))) stop("This prediction is already complete.")
  request <- readRDS(file.path(directory, "request.rds"))
  fitted <- readRDS(file.path(directory, "model.rds"))
  predictions <- fread(file.path(directory, "predictions.csv"))
  config <- request$config
  input <- readRDS(file.path(config$results_directory, "input.rds"))
  data <- soy_design(input, prediction = TRUE, canopy = request$method == "joint_canopy")
  if (!identical(names(fitted$beta), colnames(data$X))) stop("Saved coefficient columns differ from the input.")
  test <- which(data$cluster_id %in% request$test_lines)
  train <- which(data$cluster_id %in% request$training_lines)
  if (length(intersect(test, train)) || length(test) + length(train) != length(data$y))
    stop("Saved sample roles do not partition the input.")
  calibration <- setNames(vapply(config$years, function(year)
    quantile((data$y[train] - fitted$oof)[data$observations$Year[train] == year],
      request$tau, type = 1, names = FALSE), numeric(1)), config$years)
  raw <- drop(data$X[test, , drop = FALSE] %*% fitted$beta) * config$response_divisor
  calibrated <- raw + fitted$calibration[as.character(data$observations$Year[test])] * config$response_divisor
  if (!identical(predictions$cluster_id, data$cluster_id[test]) ||
      max(abs(predictions$observed - data$y[test] * config$response_divisor)) > 1e-10 ||
      max(abs(predictions$raw_prediction - raw)) > 1e-9 ||
      max(abs(predictions$prediction - calibrated)) > 1e-9 ||
      max(abs(calibration - fitted$calibration)) > 1e-12 || fitted$final_gap > 1e-8)
    stop("The saved prediction cannot be reconstructed from its fitted model.")
  loss <- (predictions$observed - predictions$prediction) * (request$tau - (predictions$observed < predictions$prediction))
  if (max(abs(loss - predictions$loss)) > 1e-10) stop("The saved pinball loss is inconsistent.")
  if (file.exists(file.path(directory, "metrics.csv"))) stop("The metric file already exists.")
  fwrite(predictions[, .(mean_loss = mean(loss), raw_mean_loss = mean(raw_loss),
    calibration = mean(covered), raw_calibration = mean(raw_covered)), by = .(Year, method, tau)],
    file.path(directory, "metrics.csv"))
  write_json(list(stage = "output_materialisation", model_and_predictions_reconstructed = TRUE,
    coefficients_refitted = FALSE, completed = as.character(Sys.time())),
    file.path(directory, "output_verification.json"), auto_unbox = TRUE, pretty = TRUE)
  write_json(list(case_id, status = "complete", method = request$method, tau = request$tau,
    training_lines = length(request$training_lines), test_lines = length(request$test_lines),
    wall_seconds = as.numeric(difftime(file.info(file.path(directory, "model.rds"))$mtime,
      as.POSIXct(request$started), units = "secs")), completed = as.character(Sys.time())),
    file.path(directory, "completion.json"), auto_unbox = TRUE, pretty = TRUE)
  backup <- file.path(config$backup_directory, "prediction", config$execution_version, case_id)
  dir.create(backup, recursive = TRUE, showWarnings = FALSE)
  for (file in list.files(directory, full.names = TRUE)) {
    target <- file.path(backup, basename(file))
    if (file.exists(target)) {
      if (file.info(file)$size != file.info(target)$size) stop("An existing backup differs from its source.")
    } else if (!file.copy(file, target)) stop("Prediction output backup failed.")
  }
  cat(case_id, "verified and completed\n")
  rm(data, input)
  gc(FALSE)
}
