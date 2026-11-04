library(data.table)
library(jsonlite)
library(yaml)
source("research/yeast_stress_application_20261104/features.R")
source("research/yeast_proteome_application_20261105/features.R")
source("research/empirical_selection_20261101/prediction_core.R")
source("research/empirical_selection_20261101/ridge_direction.R")
source("research/score_validity_20261026/load_joint.R")
source("research/joint_profile_manuscript_20261026/conic_search.R")

config <- read_yaml("config/development_ridge_offset_20261105.yml")
base <- read_yaml(config$base_config)
base$conditions <- base$assays$NaCl
directory <- config$output_root
if (dir.exists(directory)) stop("开发实验目录已经存在")
dir.create(directory, recursive = TRUE)
panel <- proteome_panel(config$input_directory)
directions <- readRDS(config$source_recipes)
rows <- fread(file.path(config$input_directory, "NaCl", "discovery.csv"))[exploration_role == "training"]
recipes <- summary <- list()
for (key in names(directions)) {
  current <- rows[strain_id %in% directions[[key]]$training_strains]
  recipe <- ep_offset_recipe(current, panel, directions[[key]], proteome_features)
  features <- ep_offset_features(current, panel, recipe, proteome_features)
  j <- recipe$penalized
  error <- max(abs(colSums(features$X[, j, drop = FALSE] * features$weight) / features$n),
    abs(colSums(features$X[, j, drop = FALSE]^2 * features$weight) / features$n - 1))
  stopifnot(error < 1e-10)
  fitted <- ep_path(features, "JOINT", 0.25, base, 2L, list(Lambda = 1, size = 2L))
  reconstruction <- ep_predict_offset(features, fitted$beta)
  direct <- features$X %*% fitted$beta * recipe$y_scale + features$ridge_prediction
  prediction_error <- max(abs(reconstruction - direct))
  stopifnot(prediction_error < 1e-12)
  recipes[[key]] <- recipe
  summary[[key]] <- data.table(panel = key, clusters = features$n, variables = ncol(features$X),
    always_included = recipe$always, residual_scale = recipe$y_scale,
    standardization_error = error, prediction_error = prediction_error,
    maximum_gap = max(vapply(fitted$diagnostics[[1L]], `[[`, numeric(1), "gap")))
}
saveRDS(recipes, file.path(directory, "recipes.rds"))
output <- "results/preflight/empirical_selection_20261101/ridge_direction"
fwrite(rbindlist(summary), file.path(output, "offset_checks.csv"))
implementation <- file.path(directory, "implementation")
dir.create(implementation)
files <- c("config/development_ridge_offset_20261105.yml", config$base_config,
           "research/empirical_selection_20261101/prediction_core.R", "research/empirical_selection_20261101/ridge_direction.R",
           "research/yeast_stress_application_20261104/features.R", "research/yeast_proteome_application_20261105/features.R",
           "research/ridge_direction_20261105/run_discovery.R")
destinations <- file.path(implementation, c("config.yml", "base_config.yml", "prediction_core.R", "ridge_direction.R",
                                           "contrasts.R", "features.R", "run_discovery.R"))
stopifnot(all(file.copy(files, destinations)))
jobs <- data.table(trait = "NaCl", method = c("JOINT_OFFSET", "POOLED_OFFSET", "JOINT_OFFSET", "JOINT_OFFSET"),
                   tau = c(0.25, 0.25, 0.5, 0.75), run_id = c("tau025", "tau025", "tau050", "tau075"))
jobs[, job_id := seq_len(.N)]
jobs[, directory := file.path(config$output_root, trait, paste0(method, "_", run_id))]
jobs[, status := "pending"]
fwrite(jobs, file.path(directory, "jobs.csv"))
print(rbindlist(summary))
