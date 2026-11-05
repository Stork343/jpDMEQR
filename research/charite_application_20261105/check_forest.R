library(data.table)
library(jsonlite)
library(yaml)
library(matrixStats)
library(grf)
source("research/charite_application_20261105/features.R")
source("research/empirical_selection_20261101/prediction_core.R")
source("research/empirical_selection_20261101/ridge_direction.R")
source("research/ridge_direction_20261105/forest.R")
config <- read_yaml("config/application_charite_hb_20261105.yml")
panel <- ch_panel(config$input_directory)
rows <- fread(file.path(config$input_directory, "training.csv"))
recipe <- readRDS(file.path(config$output_root, "recipes.rds"))$inner_1
records <- list()
for (offset in c(FALSE, TRUE)) {
  features <- function(current) if (offset) ch_offset_features(current, panel, recipe) else
    ch_features(current, panel, recipe$direction_recipe$base_recipe)
  fitting <- features(rows[inner_fold != 1L])
  held <- features(rows[inner_fold == 1L])
  settings <- list(regression_splitting = FALSE, mtry_multiplier = 1, min_node_size = 2, honesty = TRUE)
  model <- rs_forest_fit(fitting, settings, 20261105L, offset)
  prediction <- rs_forest_predict(model, held)
  native <- predict(model$model, held$X[, -1L, drop = FALSE], quantiles = unlist(config$quantiles), num.threads = 1L)$predictions
  expected <- native * held$recipe$y_scale + held$recipe$y_center
  if (offset) expected <- expected + held$offset
  error <- max(abs(prediction - expected))
  stopifnot(error < 1e-12, all(is.finite(prediction)))
  records[[length(records) + 1L]] <- data.table(offset_model = offset, patients = fitting$n, prediction_error = error,
    grf_version = model$grf_version)
}
fwrite(rbindlist(records), file.path(config$report_directory, "forest_checks.csv"))
write_json(list(status = "passed", formal_validation_used = FALSE),
  file.path(config$report_directory, "forest_checks.json"), auto_unbox = TRUE, pretty = TRUE)
