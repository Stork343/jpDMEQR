library(data.table)
library(jsonlite)
library(yaml)
library(matrixStats)
library(grf)
source("research/charite_application_20261105/features.R")
source("research/charite_application_20261105/predict.R")
source("research/empirical_selection_20261101/prediction_core.R")
source("research/ridge_direction_20261105/forest.R")
config <- read_yaml("config/application_charite_hb_20261105.yml")
root <- config$output_root
release <- read_json(file.path(root, "validation_release.json"), simplifyVector = TRUE)
stopifnot(release$status == "frozen")
columns <- c("sample_id", "strain_id", "admission_day", "age", "bmi", "charlson", "who", "condition", "role", "inner_fold")
rows <- fread(file.path(config$input_directory, "validation.csv"), select = columns)
training_ids <- fread(file.path(config$input_directory, "roles.csv"))[role == "training"]$strain_id
stopifnot(all(rows$role == "validation"), !"y" %in% names(rows), !length(intersect(rows$strain_id, training_ids)))
panel <- ch_panel(config$input_directory)
output <- file.path(root, "validation")
if (dir.exists(output)) stop("验证预测目录已经存在")
dir.create(output)
records <- list()
for (directory in release$model_directories) {
  request <- read_json(file.path(directory, "request.json"), simplifyVector = TRUE)
  if (startsWith(request$method, "XGB")) next
  fit <- readRDS(file.path(directory, "fit.rds"))
  stopifnot(setequal(fit$recipe$training_strains, training_ids))
  prediction <- ch_predict(fit, rows, panel)
  for (k in seq_along(fit$quantiles)) records[[length(records) + 1L]] <- data.table(
    method = request$method, sample_id = rows$sample_id, strain_id = rows$strain_id,
    tau = fit$quantiles[k], prediction = prediction[, k])
}
predictions <- rbindlist(records)
stopifnot(!anyDuplicated(predictions[, .(method, sample_id, tau)]))
fwrite(predictions, file.path(output, "predictions_R.csv"))
recipe <- readRDS(file.path(root, "recipes.rds"))$training
for (mode in c("raw", "offset")) {
  features <- ch_prediction_features(rows, panel,
    if (mode == "offset") recipe else recipe$direction_recipe$base_recipe, mode == "offset")
  directory <- file.path(output, "boost_features", mode)
  dir.create(directory, recursive = TRUE)
  fwrite(as.data.table(features$X[, -1L, drop = FALSE]), file.path(directory, "X.csv"))
  offset <- if (mode == "offset") features$offset else rep(0, nrow(rows))
  fwrite(data.table(sample_id = rows$sample_id, strain_id = rows$strain_id, offset = offset), file.path(directory, "metadata.csv"))
  write_json(list(y_scale = features$recipe$y_scale, y_center = features$recipe$y_center,
    columns = colnames(features$X)[-1L]), file.path(directory, "scaling.json"),
    auto_unbox = TRUE, pretty = TRUE, digits = NA)
}
write_json(list(status = "predicted", patients = uniqueN(rows$strain_id), records = nrow(rows),
  response_in_prediction_input = FALSE), file.path(output, "prediction_input_status.json"),
  auto_unbox = TRUE, pretty = TRUE)
