library(data.table)
library(jsonlite)
library(yaml)
source("research/ibdmdb_application_20261105/features.R")
source("research/empirical_selection_20261101/prediction_core.R")
source("research/empirical_selection_20261101/ridge_direction.R")
source("research/score_validity_20261026/load_joint.R")
source("research/joint_profile_manuscript_20261026/conic_search.R")
source("research/soynam_application_20261031/load_pooled_lp.R")

config <- read_yaml("config/discovery_ibdmdb_20261105.yml")
panel <- ibd_panel(config$input_directory)
rows <- fread(file.path(config$input_directory, "discovery.csv"))[outer_fold != 1L]
recipes <- readRDS(file.path(config$output_root, "fold01", "recipes.rds"))
features <- ep_offset_features(rows, panel, recipes$training, ibd_features)
press <- fread(file.path(config$output_root, "press_checks.csv"))
stopifnot(nrow(press) == 2L, all(press$press_error < 1e-7), length(unique(press$cluster_measurements)) == 2L)
checks <- rbindlist(lapply(unlist(config$quantiles), function(value) {
  fit <- ep_path(features, "JOINT", value, config, 2L, list(Lambda = 1, size = 2L))
  gap <- max(vapply(fit$diagnostics[[1L]], `[[`, numeric(1), "gap"))
  stopifnot(is.finite(gap), gap <= config$gap_tolerance, all(is.finite(fit$beta)))
  data.table(method = "JOINT", tau = value, clusters = features$n, variables = ncol(features$X), maximum_gap = gap)
}))
pooled <- ep_path(features, "POOLED", config$decision_quantile, config, 2L)
pooled_gap <- max(abs(pooled$diagnostics$primal_dual_gap))
stopifnot(is.finite(pooled_gap), pooled_gap <= config$gap_tolerance, all(is.finite(pooled$beta)))
checks <- rbind(checks, data.table(method = "POOLED", tau = config$decision_quantile,
  clusters = features$n, variables = ncol(features$X), maximum_gap = pooled_gap))
output <- "results/preflight/empirical_selection_20261101/ibdmdb_cohort"
fwrite(checks, file.path(output, "solver_checks.csv"))
fwrite(press, file.path(output, "press_checks.csv"))
write_json(list(status = "passed", formal_validation_used = FALSE),
  file.path(output, "solver_checks.json"), auto_unbox = TRUE, pretty = TRUE)
print(checks)
