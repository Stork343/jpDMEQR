library(data.table)
library(jsonlite)
library(yaml)
library(grf)
source("research/ibdmdb_application_20261105/features.R")
source("research/empirical_selection_20261101/prediction_core.R")
source("research/empirical_selection_20261101/ridge_direction.R")
source("research/ridge_direction_20261105/forest.R")
config <- read_yaml("config/development_ibdmdb_forest_20261105.yml")
panel <- ibd_panel(config$input_directory)
rows <- fread(file.path(config$input_directory, "discovery.csv"))
training <- rows[outer_fold != 1L]
evaluation <- rows[outer_fold == 1L]
recipe <- readRDS(file.path(config$source_root, "fold01", "recipes.rds"))$training
checks <- list()
for (offset in c(FALSE, TRUE)) {
  feature_function <- function(current) {
    result <- if (offset) ep_offset_features(current, panel, recipe, ibd_features) else
      ibd_features(current, panel, recipe$direction_recipe$base_recipe)
    result$sample_id <- current$sample_id
    result
  }
  fitting <- feature_function(training)
  test <- feature_function(evaluation)
  settings <- list(regression_splitting = FALSE, mtry_multiplier = 1, min_node_size = 5, honesty = TRUE)
  fit <- rs_forest_fit(fitting, settings, 20261105L, offset)
  wrapper <- rs_forest_predict(fit, test)
  native <- predict(fit$model, test$X[, -1L, drop = FALSE], quantiles = c(0.25, 0.5, 0.75), num.threads = 1L)$predictions
  expected <- native * test$recipe$y_scale + test$recipe$y_center
  if (offset) expected <- expected + test$offset
  error <- max(abs(wrapper - expected))
  stopifnot(error < 1e-12, all(wrapper[, 1L] <= wrapper[, 2L]), all(wrapper[, 2L] <= wrapper[, 3L]))
  checks[[length(checks) + 1L]] <- data.table(offset_model = offset, clusters = fitting$n,
    minimum_cluster_size = min(table(fitting$cluster_id)), maximum_cluster_size = max(table(fitting$cluster_id)),
    prediction_error = error, grf_version = fit$grf_version)
}
output <- "results/preflight/empirical_selection_20261101/ibdmdb_cohort"
fwrite(rbindlist(checks), file.path(output, "forest_checks.csv"))
write_json(list(status = "passed", formal_validation_used = FALSE),
  file.path(output, "forest_checks.json"), auto_unbox = TRUE, pretty = TRUE)
print(rbindlist(checks))
