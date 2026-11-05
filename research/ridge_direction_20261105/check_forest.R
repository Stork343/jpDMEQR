library(data.table)
library(jsonlite)
library(yaml)
library(grf)
source("research/yeast_stress_application_20261104/features.R")
source("research/yeast_proteome_application_20261105/features.R")
source("research/empirical_selection_20261101/ridge_direction.R")
source("research/ridge_direction_20261105/forest.R")

config <- read_yaml("config/development_ridge_offset_series_20261105.yml")
panel <- proteome_panel(config$input_directory)
recipes <- readRDS(file.path(config$output_root, "NaCl", "recipes.rds"))
rows <- fread(file.path(config$input_directory, "NaCl", "discovery.csv"))[exploration_role == "training"]
ids <- sort(unique(rows$strain_id))[unique(as.integer(seq(1, uniqueN(rows$strain_id), length.out = 60L)))]
sample <- rows[strain_id %in% ids]
features <- ep_offset_features(sample, panel, recipes$training, proteome_features)
features$sample_id <- sample$sample_id
keep <- unique(c(1L:12L, as.integer(seq(13, ncol(features$X), length.out = 30L))))
features$X <- features$X[, keep, drop = FALSE]
settings <- list(regression_splitting = FALSE, mtry_multiplier = 1, min_node_size = 5L, honesty = TRUE)
first <- rs_forest_fit(features, settings, 20261105L, TRUE, compute_oob = TRUE)
second <- rs_forest_fit(features, settings, 20261105L, TRUE, compute_oob = FALSE)
first_prediction <- rs_forest_predict(first, features)
second_prediction <- rs_forest_predict(second, features)
model_error <- max(abs(first_prediction - second_prediction))
native <- predict(first$model, features$X[, -1L, drop = FALSE], quantiles = c(0.25, 0.5, 0.75), num.threads = 1L)$predictions
expected <- native * features$recipe$y_scale + features$recipe$y_center + features$offset
mapping_error <- max(abs(first_prediction - expected))
oob <- rs_forest_predict(first, features, oob = TRUE)
stopifnot(model_error < 1e-12, mapping_error < 1e-12, all(is.finite(oob)))
output <- "results/preflight/empirical_selection_20261101/ridge_direction"
write_json(list(status = "passed", actual_training_clusters = features$n,
  trees = 2000L, oob_option_prediction_error = model_error, native_prediction_mapping_error = mapping_error,
  package_version = as.character(packageVersion("grf")), formal_validation_used = FALSE),
  file.path(output, "forest_checks.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
cat("Quantile forest 输入与预测换算检查通过\n")
