library(data.table)
library(jsonlite)
library(yaml)
source("research/yeast_stress_application_20261104/features.R")
source("research/yeast_proteome_application_20261105/features.R")
source("research/empirical_selection_20261101/prediction_core.R")
source("research/score_validity_20261026/load_joint.R")
source("research/joint_profile_manuscript_20261026/conic_search.R")

config <- read_yaml("config/discovery_yeast_proteome_20261105.yml")
source_check <- read_json("results/preflight/empirical_selection_20261101/yeast_proteome/source_checks.json")
stopifnot(source_check$status == "passed")
panel <- proteome_panel(config$input_directory)
records <- list()
for (assay_name in c("Anisomycin", "NaCl")) {
  rows <- fread(file.path(config$input_directory, assay_name, "discovery.csv"))[exploration_role == "training"]
  current <- config
  current$conditions <- config$assays[[assay_name]]
  recipe <- proteome_recipe(rows, panel, current)
  features <- proteome_features(rows, panel, recipe)
  d <- length(recipe$low_columns)
  p <- length(recipe$genes)
  P <- features$X[, d + seq_len(p), drop = FALSE]
  center_error <- max(abs(colMeans(P)))
  scale_error <- max(abs(colMeans(P^2) - 1))
  condition <- recipe$contrasts[rows$condition, , drop = FALSE]
  error <- 0
  for (j in seq_len(ncol(condition))) error <- max(error,
    abs(features$X[, d + j * p + seq_len(p)] - P * condition[, j]))
  stopifnot(center_error < 1e-10, scale_error < 1e-10, error < 1e-12,
            max(abs(rowsum(features$weight, features$cluster_id) - 1)) < 1e-12)
  fit <- ep_path(features, "JOINT", 0.25, config, 2L, list(Lambda = 1, size = 2L))
  gaps <- vapply(fit$diagnostics[[1L]], `[[`, numeric(1), "gap")
  stopifnot(max(gaps) <= config$gap_tolerance)
  records[[assay_name]] <- data.table(assay = assay_name, strains = features$n, proteins = p,
    variables = ncol(features$X), adjustment_variables = d, center_error = center_error,
    scale_error = scale_error, interaction_error = error, maximum_gap = max(gaps))
}
output <- "results/preflight/empirical_selection_20261101/yeast_proteome"
fwrite(rbindlist(records), file.path(output, "input_checks.csv"))
write_json(list(status = "passed", formal_validation_used = FALSE, manuscript_adoption = FALSE),
           file.path(output, "input_checks.json"), auto_unbox = TRUE, pretty = TRUE)
print(rbindlist(records))
