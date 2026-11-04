library(data.table)
library(jsonlite)
library(yaml)
source("research/yeast_stress_application_20261104/features.R")
source("research/yeast_proteome_application_20261105/features.R")
source("research/empirical_selection_20261101/prediction_core.R")
source("research/empirical_selection_20261101/ridge_direction.R")

config <- read_yaml("config/development_ridge_direction_20261105.yml")
base <- read_yaml(config$base_config)
base$conditions <- base$assays[[config$assay]]
checks <- read_json("results/preflight/empirical_selection_20261101/ridge_direction/checks.json", simplifyVector = TRUE)
stopifnot(checks$status == "passed")
directory <- config$output_root
if (dir.exists(directory)) stop("开发实验目录已经存在")
dir.create(directory, recursive = TRUE)
panel <- proteome_panel(base$input_directory)
rows <- fread(file.path(base$input_directory, config$assay, "discovery.csv"))[exploration_role == "training"]
recipes <- list()
summary <- list()
for (key in c("training", "inner_1", "inner_2", "inner_3")) {
  current <- if (key == "training") rows else rows[inner_fold != as.integer(sub("inner_", "", key))]
  started <- proc.time()[["elapsed"]]
  recipe <- ep_direction_recipe(current, panel, base, proteome_recipe, proteome_features)
  features <- ep_direction_features(current, panel, recipe, proteome_features)
  directions <- tail(seq_len(ncol(features$X)), length(recipe$direction_scale))
  error <- max(abs(colSums(features$X[, directions, drop = FALSE]^2 * features$weight) / features$n - 1))
  stopifnot(error < 1e-10)
  recipes[[key]] <- recipe
  summary[[key]] <- data.table(panel = key, clusters = features$n,
    ridge_lambda = recipe$background$lambda[recipe$background$choice],
    minimum_delete_margin = recipe$background$minimum_margin,
    system_error = recipe$background$system_error, direction_scale_error = error,
    runtime_sec = proc.time()[["elapsed"]] - started)
  cat(key, "completed\n")
  flush.console()
}
saveRDS(recipes, file.path(directory, "recipes.rds"))
fwrite(rbindlist(summary), file.path(directory, "recipe_summary.csv"))
output <- "results/preflight/empirical_selection_20261101/ridge_direction"
fwrite(rbindlist(summary), file.path(output, "recipe_summary.csv"))
implementation <- file.path(directory, "implementation")
dir.create(implementation)
files <- c("config/development_ridge_direction_20261105.yml", config$base_config,
           "research/empirical_selection_20261101/prediction_core.R", "research/empirical_selection_20261101/ridge_direction.R",
           "research/yeast_stress_application_20261104/features.R", "research/yeast_proteome_application_20261105/features.R",
           "research/ridge_direction_20261105/run_discovery.R")
destinations <- file.path(implementation, c("config.yml", "base_config.yml", "prediction_core.R", "ridge_direction.R",
                                           "contrasts.R", "features.R", "run_discovery.R"))
stopifnot(all(file.copy(files, destinations)))
jobs <- data.table(trait = "NaCl", method = c("JOINT_DIRECTION", "POOLED_DIRECTION", "PRESS_RIDGE", "JOINT_DIRECTION", "JOINT_DIRECTION"),
                   tau = c(0.25, 0.25, NA_real_, 0.5, 0.75), run_id = c("tau025", "tau025", "cv", "tau050", "tau075"))
jobs[, job_id := seq_len(.N)]
jobs[, directory := file.path(config$output_root, trait, paste0(method, "_", run_id))]
jobs[, status := "pending"]
fwrite(jobs, file.path(directory, "jobs.csv"))
