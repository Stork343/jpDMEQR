library(data.table)
library(jsonlite)
library(yaml)
library(BEDMatrix)
source("research/arabidopsis_application_20261105/features.R")
source("research/empirical_selection_20261101/prediction_core.R")
source("research/score_validity_20261026/load_joint.R")
source("research/joint_profile_manuscript_20261026/conic_search.R")

config <- read_yaml("config/discovery_arabidopsis_20261105.yml")
panel <- arab_panel(config$input_directory)
rows <- fread(file.path(config$input_directory, "discovery.csv"), colClasses = c(accession_id = "character"))[exploration_role == "training"]
recipe <- arab_recipe(rows, panel, config)
ids <- recipe$training_accessions[unique(as.integer(seq(1, length(recipe$training_accessions), length.out = 40L)))]
features <- arab_features(rows[accession_id %in% ids], panel, recipe)
d <- length(recipe$low_columns)
p <- length(recipe$markers)
selected <- unique(as.integer(seq(1, p, length.out = 64L)))
keep <- c(seq_len(d), d + selected, d + p + selected)
features$X <- features$X[, keep, drop = FALSE]
features$penalty_factor <- features$penalty_factor[keep]
results <- rbindlist(lapply(unlist(config$quantiles), function(value) {
  fitted <- ep_path(features, "JOINT", value, config, 2L, list(Lambda = 1, size = 2L))
  gap <- max(vapply(fitted$diagnostics[[1L]], `[[`, numeric(1), "gap"))
  stopifnot(is.finite(gap), gap <= config$gap_tolerance, all(is.finite(fitted$beta)))
  data.table(tau = value, actual_training_accessions = features$n, variables = ncol(features$X), maximum_gap = gap)
}))
output <- "results/preflight/empirical_selection_20261101/arabidopsis"
fwrite(results, file.path(output, "solver_checks.csv"))
write_json(list(status = "passed", formal_validation_used = FALSE),
           file.path(output, "solver_checks.json"), auto_unbox = TRUE, pretty = TRUE)
print(results)
