library(data.table)
library(jsonlite)
library(yaml)
library(matrixStats)
config <- read_yaml("config/development_drugmatrix_offset_20261105.yml")
implementation <- file.path(config$output_root, "implementation")
for (file in c("dm_features.R", "offset_features.R", "prediction_core.R", "ridge_direction.R")) source(file.path(implementation, file))
source("research/empirical_selection_20261101/report_prediction.R")
queue <- fread(file.path(config$output_root, "progress.csv"))
stopifnot(nrow(queue) == 7L, all(queue$status == "completed"))
panel <- dmo_panel(config)
offset_features <- function(rows, panel, recipe) ep_offset_features(rows, panel, recipe, dmo_features)
feature_functions <- list(JOINT_OFFSET = offset_features, POOLED_OFFSET = offset_features,
  JOINT_DIRECT = dmo_features, ENET = dmo_features, PRESS_RIDGE = dmo_features)
prediction_functions <- list(JOINT_OFFSET = ep_predict_offset, POOLED_OFFSET = ep_predict_offset,
  JOINT_DIRECT = ep_predict, ENET = ep_predict, PRESS_RIDGE = ep_predict)
output <- "results/preflight/empirical_selection_20261101/drugmatrix_offset"
comparison <- ep_report_prediction(config, queue, panel, dmo_features, output, target_method = "JOINT_OFFSET",
  stage = "method_development", feature_functions = feature_functions, prediction_functions = prediction_functions,
  cluster_count_column = "compounds")
print(comparison)
