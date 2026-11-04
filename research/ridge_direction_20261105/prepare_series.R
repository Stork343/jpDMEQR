library(data.table)
library(jsonlite)
library(yaml)
source("research/yeast_stress_application_20261104/features.R")
source("research/yeast_proteome_application_20261105/features.R")
source("research/empirical_selection_20261101/prediction_core.R")
source("research/empirical_selection_20261101/ridge_direction.R")

config <- read_yaml("config/development_ridge_offset_series_20261105.yml")
base <- read_yaml(config$base_config)
check <- read_json("results/preflight/empirical_selection_20261101/ridge_direction/three_condition_checks.json")
stopifnot(check$status == "passed")
root <- config$output_root
if (dir.exists(root)) stop("开发实验目录已经存在")
dir.create(root, recursive = TRUE)
panel <- proteome_panel(config$input_directory)
summary <- jobs <- list()
for (assay_name in unlist(config$assays)) {
  directory <- file.path(root, assay_name)
  dir.create(directory)
  base$conditions <- base$assays[[assay_name]]
  rows <- fread(file.path(config$input_directory, assay_name, "discovery.csv"))[exploration_role == "training"]
  recipes <- list()
  for (key in c("training", "inner_1", "inner_2", "inner_3")) {
    current <- if (key == "training") rows else rows[inner_fold != as.integer(sub("inner_", "", key))]
    started <- proc.time()[["elapsed"]]
    direction <- ep_direction_recipe(current, panel, base, proteome_recipe, proteome_features)
    recipe <- ep_offset_recipe(current, panel, direction, proteome_features)
    features <- ep_offset_features(current, panel, recipe, proteome_features)
    j <- recipe$penalized
    error <- max(abs(colSums(features$X[, j, drop = FALSE] * features$weight) / features$n),
      abs(colSums(features$X[, j, drop = FALSE]^2 * features$weight) / features$n - 1))
    stopifnot(error < 1e-10)
    recipes[[key]] <- recipe
    summary[[length(summary) + 1L]] <- data.table(assay = assay_name, panel = key, clusters = features$n,
      variables = ncol(features$X), residual_scale = recipe$y_scale, standardization_error = error,
      ridge_lambda = direction$background$lambda[direction$background$choice],
      background_and_transform_sec = proc.time()[["elapsed"]] - started)
  }
  saveRDS(recipes, file.path(directory, "recipes.rds"))
  jobs[[assay_name]] <- data.table(trait = assay_name,
    method = c("JOINT_OFFSET", "PRESS_RIDGE", "POOLED_OFFSET", "JOINT_OFFSET", "JOINT_OFFSET"),
    tau = c(0.25, NA_real_, 0.25, 0.5, 0.75), run_id = c("tau025", "cv", "tau025", "tau050", "tau075"))
  cat(assay_name, "prepared\n")
  flush.console()
}
fwrite(rbindlist(summary), file.path(root, "recipe_summary.csv"))
implementation <- file.path(root, "implementation")
dir.create(implementation)
files <- c("config/development_ridge_offset_series_20261105.yml", config$base_config,
           "research/empirical_selection_20261101/prediction_core.R", "research/empirical_selection_20261101/ridge_direction.R",
           "research/yeast_stress_application_20261104/features.R", "research/yeast_proteome_application_20261105/features.R",
           "research/ridge_direction_20261105/run_discovery.R")
destinations <- file.path(implementation, c("config.yml", "base_config.yml", "prediction_core.R", "ridge_direction.R",
                                           "contrasts.R", "features.R", "run_discovery.R"))
stopifnot(all(file.copy(files, destinations)))
queue <- rbindlist(jobs)
queue[, job_id := seq_len(.N)]
queue[, directory := file.path(root, trait, paste0(method, "_", run_id))]
queue[, status := "pending"]
fwrite(queue, file.path(root, "jobs.csv"))
output <- "results/preflight/empirical_selection_20261101/ridge_direction"
fwrite(rbindlist(summary), file.path(output, "series_recipe_summary.csv"))
