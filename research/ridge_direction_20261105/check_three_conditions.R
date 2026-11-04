library(data.table)
library(jsonlite)
library(yaml)
source("research/yeast_stress_application_20261104/features.R")
source("research/yeast_proteome_application_20261105/features.R")
source("research/empirical_selection_20261101/prediction_core.R")
source("research/empirical_selection_20261101/ridge_direction.R")
config <- read_yaml("config/discovery_yeast_proteome_20261105.yml")
config$conditions <- config$assays$Anisomycin
panel <- proteome_panel(config$input_directory)
rows <- fread(file.path(config$input_directory, "Anisomycin", "discovery.csv"))[exploration_role == "training"]
recipe <- proteome_recipe(rows, panel, config)
ids <- sort(unique(rows$strain_id))[unique(as.integer(seq(1, uniqueN(rows$strain_id), length.out = 60L)))]
sample <- rows[strain_id %in% ids]
features <- proteome_features(sample, panel, recipe)
low <- which(features$penalty_factor == 0)
high <- which(features$penalty_factor == 1)
keep <- c(low, high[unique(as.integer(seq(1, length(high), length.out = 30L)))])
features$X <- features$X[, keep, drop = FALSE]
features$penalty_factor <- features$penalty_factor[keep]
config$elastic_net_path_length <- 5L
fit <- ep_ridge_press(features, config)
checks <- list()
for (id in head(ids, 3L)) {
  held <- which(features$cluster_id == id)
  retained <- setdiff(seq_len(nrow(features$X)), held)
  reference <- features
  reference$X <- features$X[retained, , drop = FALSE]
  reference$y <- features$y[retained]
  reference$weight <- features$weight[retained]
  reference$cluster_id <- features$cluster_id[retained]
  reference$n <- features$n - 1L
  model <- ep_mean_path(reference, 0, config, lambda = fit$lambda * features$n / reference$n)
  prediction <- as.matrix(features$X[held, , drop = FALSE] %*% model$beta)
  error <- max(abs(prediction - fit$leave_out_path[held, , drop = FALSE]))
  stopifnot(error < 1e-8)
  checks[[id]] <- data.table(cluster = id, conditions = length(held), prediction_error = error)
}
output <- "results/preflight/empirical_selection_20261101/ridge_direction"
fwrite(rbindlist(checks), file.path(output, "three_condition_checks.csv"))
write_json(list(status = "passed", formal_validation_used = FALSE),
           file.path(output, "three_condition_checks.json"), auto_unbox = TRUE, pretty = TRUE)
print(rbindlist(checks))
