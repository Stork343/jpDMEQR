library(data.table)
library(jsonlite)
library(yaml)
library(glmnet)
source("research/yeast_stress_application_20261104/features.R")
source("research/yeast_proteome_application_20261105/features.R")
source("research/empirical_selection_20261101/prediction_core.R")
source("research/empirical_selection_20261101/ridge_direction.R")

config <- read_yaml("config/discovery_yeast_proteome_20261105.yml")
config$conditions <- config$assays$NaCl
panel <- proteome_panel(config$input_directory)
rows <- fread(file.path(config$input_directory, "NaCl", "discovery.csv"))[exploration_role == "training"]
recipe <- proteome_recipe(rows, panel, config)
sample_ids <- unique(rows$strain_id)[unique(as.integer(seq(1, uniqueN(rows$strain_id), length.out = 60L)))]
sample <- rows[strain_id %in% sample_ids]
features <- proteome_features(sample, panel, recipe)
low <- which(features$penalty_factor == 0)
high <- which(features$penalty_factor == 1)
keep <- c(low, high[unique(as.integer(seq(1, length(high), length.out = 30L)))])
features$X <- features$X[, keep, drop = FALSE]
features$penalty_factor <- features$penalty_factor[keep]
config$elastic_net_path_length <- 5L
fit <- ep_ridge_press(features, config)
reference <- ep_mean_path(features, 0, config, lambda = fit$lambda)
coefficient_error <- max(abs(fit$beta - as.numeric(reference$beta[, fit$choice])))
stopifnot(coefficient_error < 1e-8)
checks <- list()
for (id in head(sample_ids, 3L)) {
  held <- which(features$cluster_id == id)
  retained <- setdiff(seq_len(nrow(features$X)), held)
  reference_features <- features
  reference_features$X <- features$X[retained, , drop = FALSE]
  reference_features$y <- features$y[retained]
  reference_features$weight <- features$weight[retained]
  reference_features$cluster_id <- features$cluster_id[retained]
  reference_features$n <- features$n - 1L
  reference_fit <- ep_mean_path(reference_features, 0, config,
    lambda = fit$lambda * features$n / reference_features$n)
  predicted <- as.matrix(features$X[held, , drop = FALSE] %*% reference_fit$beta)
  error <- max(abs(predicted - fit$leave_out_path[held, , drop = FALSE]))
  stopifnot(error < 1e-8)
  checks[[id]] <- data.table(training_cluster = id, actual_records = length(held), maximum_prediction_error = error)
}
output <- "results/preflight/empirical_selection_20261101/ridge_direction"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
fwrite(rbindlist(checks), file.path(output, "leave_cluster_out_checks.csv"))
write_json(list(status = "passed", coefficient_error = coefficient_error, system_error = fit$system_error,
                minimum_leave_out_margin = fit$minimum_margin, actual_training_clusters = features$n,
                formal_validation_used = FALSE), file.path(output, "checks.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
print(rbindlist(checks))
