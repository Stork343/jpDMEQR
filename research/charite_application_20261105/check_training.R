library(data.table)
library(jsonlite)
library(yaml)
library(matrixStats)
source("research/charite_application_20261105/features.R")
source("research/charite_application_20261105/train_only.R")
source("research/empirical_selection_20261101/prediction_core.R")
source("research/empirical_selection_20261101/ridge_direction.R")
source("research/score_validity_20261026/load_joint.R")
source("research/joint_profile_manuscript_20261026/conic_search.R")
source("research/soynam_application_20261031/load_pooled_lp.R")
config <- read_yaml("config/application_charite_hb_20261105.yml")
panel <- ch_panel(config$input_directory)
rows <- fread(file.path(config$input_directory, "training.csv"))
recipes <- readRDS(file.path(config$output_root, "recipes.rds"))
features <- ch_offset_features(rows, panel, recipes$training)
checks <- rbindlist(lapply(unlist(config$quantiles), function(value) {
  fit <- ep_path(features, "JOINT", value, config, 2L, list(Lambda = 1, size = 2L))
  gap <- max(vapply(fit$diagnostics[[1L]], `[[`, numeric(1), "gap"))
  stopifnot(gap <= config$gap_tolerance, all(is.finite(fit$beta)))
  data.table(check = "joint_qp", tau = value, maximum_error = gap)
}))
base <- ch_features(rows, panel, recipes$training$direction_recipe$base_recipe)
background <- recipes$training$direction_recipe$background
sizes <- table(base$cluster_id)
for (id in unique(c(names(sizes)[which.min(sizes)], names(sizes)[which.max(sizes)]))) {
  keep <- base$cluster_id != id
  reference <- base
  for (field in c("X", "Z")) reference[[field]] <- base[[field]][keep, , drop = FALSE]
  for (field in c("y", "raw_y", "cluster_id", "weight")) reference[[field]] <- base[[field]][keep]
  reference$n <- base$n - 1L
  path <- ep_mean_path(reference, 0, config, background$lambda * base$n / reference$n)
  prediction <- drop(base$X[!keep, , drop = FALSE] %*% path$beta[, background$choice])
  error <- max(abs(prediction - background$leave_out_prediction[rows$sample_id[!keep]]))
  stopifnot(error < 1e-7)
  checks <- rbind(checks, data.table(check = paste0("ridge_delete_", sum(!keep)), tau = NA_real_, maximum_error = error))
}
fitting <- rows[inner_fold != 1L]
evaluation <- rows[inner_fold == 1L]
fitting[, inner_fold := rep(1:3, length.out = uniqueN(strain_id))[match(strain_id, sort(unique(strain_id)))]]
progress <- function(stage, details = list()) invisible(NULL)
root <- file.path(config$output_root, "checks")
dir.create(root)
for (method in c("JOINT", "POOLED", "ENET")) {
  directories <- file.path(root, paste0(method, c("_training", "_reference")))
  stopifnot(all(vapply(directories, dir.create, logical(1))))
  mode <- if (method == "ENET") "oof" else "training"
  learned <- ch_train_prediction(fitting, panel, method, 0.25, config, progress, directories[1L],
    ch_recipe, ch_features, ep_predict, mode)
  reference <- ep_fit_prediction(fitting, evaluation, panel, method, 0.25, config, progress, directories[2L],
    ch_recipe, ch_features, "strain_id", "condition", ep_predict, mode)
  error <- max(abs(learned$beta - reference$beta), abs(unlist(learned$calibration) - unlist(reference$calibration)))
  stopifnot(error < 1e-10, identical(learned$choice, reference$choice))
  checks <- rbind(checks, data.table(check = paste0("training_interface_", method), tau = 0.25, maximum_error = error))
}
output <- config$report_directory
fwrite(checks, file.path(output, "training_checks.csv"))
write_json(list(status = "passed", formal_validation_used = FALSE), file.path(output, "training_checks.json"),
  auto_unbox = TRUE, pretty = TRUE)
print(checks)
