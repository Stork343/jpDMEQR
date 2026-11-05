library(data.table)
library(jsonlite)
library(yaml)
library(matrixStats)
source("research/drugmatrix_application_20261101/prediction.R")
source("research/drugmatrix_application_20261101/offset_features.R")
source("research/empirical_selection_20261101/prediction_core.R")
source("research/empirical_selection_20261101/ridge_direction.R")
config_path <- "config/development_drugmatrix_offset_20261105.yml"
config <- read_yaml(config_path)
root <- config$output_root
if (dir.exists(root)) stop("运行目录已经存在")
dir.create(root, recursive = TRUE)
rows <- readRDS(config$input_growth)
stopifnot(all(rows$role == "discovery"), all(is.finite(rows$y)))
rows[, `:=`(strain_id = compound_id, condition = "ALT")]
training <- rows[exploration_role == "training"]
evaluation <- rows[exploration_role == "evaluation"]
stopifnot(!length(intersect(training$strain_id, evaluation$strain_id)))
derived <- file.path(config$input_directory, "ALT")
if (dir.exists(derived)) stop("开发数据目录已经存在")
dir.create(derived, recursive = TRUE)
fwrite(rows, file.path(derived, "discovery.csv"))
panel <- dmo_panel(config)
recipes <- list()
summary <- checks <- list()
for (key in c("training", "inner_1", "inner_2", "inner_3")) {
  current <- if (key == "training") training else training[inner_fold != as.integer(sub("inner_", "", key))]
  started <- proc.time()[["elapsed"]]
  direction <- ep_direction_recipe(current, panel, config, dmo_recipe, dmo_features)
  recipe <- ep_offset_recipe(current, panel, direction, dmo_features)
  features <- ep_offset_features(current, panel, recipe, dmo_features)
  j <- recipe$penalized
  error <- max(abs(colSums(features$X[, j, drop = FALSE] * features$weight) / features$n),
    abs(colSums(features$X[, j, drop = FALSE]^2 * features$weight) / features$n - 1))
  stopifnot(error < 1e-10)
  recipes[[key]] <- recipe
  summary[[length(summary) + 1L]] <- data.table(panel = key, compounds = features$n,
    animals = nrow(features$X), variables = ncol(features$X), standardization_error = error,
    background_and_transform_sec = proc.time()[["elapsed"]] - started)
  if (key == "training") {
    base <- dmo_features(current, panel, direction$base_recipe)
    sizes <- table(base$cluster_id)
    for (id in unique(c(names(sizes)[which.min(sizes)], names(sizes)[which.max(sizes)]))) {
      keep <- base$cluster_id != id
      reference <- base
      for (field in c("X", "Z")) reference[[field]] <- base[[field]][keep, , drop = FALSE]
      for (field in c("y", "raw_y", "cluster_id", "weight")) reference[[field]] <- base[[field]][keep]
      reference$n <- base$n - 1L
      path <- ep_mean_path(reference, 0, config, direction$background$lambda * base$n / reference$n)
      prediction <- drop(base$X[!keep, , drop = FALSE] %*% path$beta[, direction$background$choice])
      press_error <- max(abs(prediction - direction$background$leave_out_prediction[current$sample_id[!keep]]))
      stopifnot(press_error < 1e-7)
      checks[[length(checks) + 1L]] <- data.table(animals = sum(!keep), press_error)
      rm(reference, path)
      gc()
    }
    rm(base)
  }
  rm(features, direction)
  gc()
  cat(key, "prepared\n")
  flush.console()
}
dir.create(file.path(root, "ALT"))
saveRDS(recipes, file.path(root, "ALT", "recipes.rds"))
fwrite(rbindlist(summary), file.path(root, "recipe_summary.csv"))
fwrite(rbindlist(checks), file.path(root, "press_checks.csv"))
implementation <- file.path(root, "implementation")
dir.create(implementation)
paths <- c(config_path, "research/drugmatrix_application_20261101/prediction.R",
  "research/drugmatrix_application_20261101/offset_features.R", "research/empirical_selection_20261101/prediction_core.R",
  "research/empirical_selection_20261101/ridge_direction.R", "research/empirical_selection_20261101/run_prediction.R",
  "research/drugmatrix_application_20261101/run_offset.R")
stopifnot(all(file.copy(paths, file.path(implementation, c("config.yml", "dm_features.R", "offset_features.R",
  "prediction_core.R", "ridge_direction.R", "run_prediction.R", "run_discovery.R")))))
queue <- data.table(trait = "ALT",
  method = c("JOINT_OFFSET", "PRESS_RIDGE", "ENET", "POOLED_OFFSET", "JOINT_DIRECT", "JOINT_OFFSET", "JOINT_OFFSET"),
  tau = c(0.75, NA_real_, NA_real_, 0.75, 0.75, 0.25, 0.5),
  run_id = c("tau075", "cv", "cv", "tau075", "tau075", "tau025", "tau050"))
queue[, job_id := seq_len(.N)]
queue[, directory := file.path(root, trait, paste0(method, "_", run_id))]
queue[, status := "pending"]
fwrite(queue, file.path(root, "jobs.csv"))
