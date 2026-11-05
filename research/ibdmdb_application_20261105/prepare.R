library(data.table)
library(jsonlite)
library(yaml)
source("research/ibdmdb_application_20261105/features.R")
source("research/empirical_selection_20261101/prediction_core.R")
source("research/empirical_selection_20261101/ridge_direction.R")

config_path <- "config/discovery_ibdmdb_20261105.yml"
config <- read_yaml(config_path)
audit <- read_json("results/preflight/empirical_selection_20261101/ibdmdb_cohort/status.json", simplifyVector = TRUE)
stopifnot(audit$status == "prepared", !audit$outcome_associations_fitted)
root <- config$output_root
if (dir.exists(root)) stop("运行目录已经存在")
dir.create(root, recursive = TRUE)
panel <- ibd_panel(config$input_directory)
rows <- fread(file.path(config$input_directory, "discovery.csv"))
stopifnot(all(rows$role == "discovery"), identical(sort(unique(rows$outer_fold)), 1:5))
jobs <- checks <- preparation <- list()
for (fold in 1:5) {
  directory <- file.path(root, sprintf("fold%02d", fold))
  dir.create(directory)
  training <- copy(rows[outer_fold != fold])
  participants <- unique(training[, .(strain_id, site)])
  setorder(participants, site, strain_id)
  set.seed(20261105L + fold)
  participants[, inner_fold := sample(rep(1:3, length.out = .N)), by = site]
  training[, inner_fold := participants$inner_fold[match(strain_id, participants$strain_id)]]
  fwrite(training[, .(strain_id, sample_id, inner_fold)], file.path(directory, "training_roles.csv"))
  recipes <- list()
  for (key in c("training", "inner_1", "inner_2", "inner_3")) {
    current <- if (key == "training") training else training[inner_fold != as.integer(sub("inner_", "", key))]
    started <- proc.time()[["elapsed"]]
    direction <- ep_direction_recipe(current, panel, config, ibd_recipe, ibd_features)
    recipe <- ep_offset_recipe(current, panel, direction, ibd_features)
    features <- ep_offset_features(current, panel, recipe, ibd_features)
    j <- recipe$penalized
    error <- max(abs(colSums(features$X[, j, drop = FALSE] * features$weight) / features$n),
      abs(colSums(features$X[, j, drop = FALSE]^2 * features$weight) / features$n - 1))
    stopifnot(error < 1e-10)
    recipes[[key]] <- recipe
    preparation[[length(preparation) + 1L]] <- data.table(fold, panel = key, clusters = features$n,
      variables = ncol(features$X), standardization_error = error,
      elapsed_sec = proc.time()[["elapsed"]] - started)
    if (fold == 1L && key == "training") {
      base <- ibd_features(current, panel, direction$base_recipe)
      sizes <- table(base$cluster_id)
      for (id in unique(c(names(sizes)[which.min(sizes)], names(sizes)[which.max(sizes)]))) {
        included <- base$cluster_id != id
        reference <- base
        for (field in c("X", "Z")) reference[[field]] <- base[[field]][included, , drop = FALSE]
        for (field in c("y", "raw_y", "cluster_id", "weight")) reference[[field]] <- base[[field]][included]
        reference$n <- base$n - 1L
        path <- ep_mean_path(reference, 0, config, direction$background$lambda * base$n / reference$n)
        prediction <- drop(base$X[!included, , drop = FALSE] %*% path$beta[, direction$background$choice])
        press_error <- max(abs(prediction - direction$background$leave_out_prediction[current$sample_id[!included]]))
        stopifnot(press_error < 1e-7)
        checks[[length(checks) + 1L]] <- data.table(fold, cluster_measurements = sum(!included), press_error)
      }
    }
  }
  saveRDS(recipes, file.path(directory, "recipes.rds"))
  jobs[[fold]] <- data.table(trait = sprintf("fold%02d", fold),
    method = c("JOINT_OFFSET", "PRESS_RIDGE", "ENET", "POOLED_OFFSET", "JOINT_DIRECT", "JOINT_OFFSET", "JOINT_OFFSET"),
    tau = c(0.75, NA_real_, NA_real_, 0.75, 0.75, 0.25, 0.5),
    run_id = c("tau075", "cv", "cv", "tau075", "tau075", "tau025", "tau050"))
  cat("Prepared", fold, "\n")
  flush.console()
}
implementation <- file.path(root, "implementation")
dir.create(implementation)
paths <- c(config_path, "research/empirical_selection_20261101/prediction_core.R",
  "research/empirical_selection_20261101/ridge_direction.R", "research/ibdmdb_application_20261105/features.R",
  "research/ibdmdb_application_20261105/run_case.R")
stopifnot(all(file.copy(paths, file.path(implementation,
  c("config.yml", "prediction_core.R", "ridge_direction.R", "features.R", "run_discovery.R")))))
queue <- rbindlist(jobs)
queue[, job_id := seq_len(.N)]
queue[, directory := file.path(root, trait, paste0(method, "_", run_id))]
queue[, status := "pending"]
fwrite(queue, file.path(root, "jobs.csv"))
fwrite(rbindlist(checks), file.path(root, "press_checks.csv"))
stopifnot(length(checks) == 2L)
fwrite(rbindlist(preparation), file.path(root, "recipe_summary.csv"))
write_json(list(status = "prepared", formal_validation_used = FALSE, manuscript_adoption = FALSE),
  file.path(root, "status.json"), auto_unbox = TRUE, pretty = TRUE)
