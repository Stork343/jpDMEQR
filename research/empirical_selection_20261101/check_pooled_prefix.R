library(data.table)
library(jsonlite)
library(yaml)
library(matrixStats)
library(filelock)
library(Rcpp)
library(RcppEigen)
source("research/soynam_application_20261031/load_pooled_lp.R")
source("research/empirical_selection_20261101/prediction_core.R")
source("research/empirical_selection_20261101/ridge_direction.R")
source("research/ibdmdb_application_20261105/features.R")
source("research/drugmatrix_application_20261101/prediction.R")
source("research/drugmatrix_application_20261101/offset_features.R")
records <- list()

for (dataset in c("IBDMDB", "DrugMatrix")) {
  if (dataset == "IBDMDB") {
    config <- read_yaml("config/discovery_ibdmdb_20261105.yml")
    panel <- ibd_panel(config$input_directory)
    rows <- fread(file.path(config$input_directory, "discovery.csv"))[outer_fold != 1L]
    roles <- fread(file.path(config$output_root, "fold01", "training_roles.csv"))
    stopifnot(identical(rows$sample_id, roles$sample_id))
    rows[, inner_fold := roles$inner_fold]
    fixture <- readRDS(file.path(config$output_root, "fold01", "POOLED_OFFSET_tau075", "inner_1.rds"))
    features <- ep_offset_features(rows[inner_fold != 1L], panel, fixture$recipe, ibd_features)
    endpoints <- c(1L, 5L, 17L, 40L)
  } else {
    config <- read_yaml("config/development_drugmatrix_offset_20261105.yml")
    panel <- dmo_panel(config)
    rows <- fread(file.path(config$input_directory, "ALT", "discovery.csv"))[exploration_role == "training"]
    fixture <- readRDS(file.path(config$output_root, "ALT", "POOLED_OFFSET_tau075", "inner_1.rds"))
    features <- ep_offset_features(rows[inner_fold != 1L], panel, fixture$recipe, dmo_features)
    endpoints <- c(1L, 5L)
  }
  row_weight <- nrow(features$X) * features$weight / features$n
  for (endpoint in endpoints) {
    started <- proc.time()[["elapsed"]]
    fit <- soy_pooled_lp_cpp(features$X * row_weight, features$y * row_weight,
      features$penalty_factor, 0.75, 40L, 0.01, numeric(), FALSE, maximum_solutions = endpoint)
    elapsed <- proc.time()[["elapsed"]] - started
    reference <- fixture$path$beta[, seq_len(endpoint), drop = FALSE]
    coefficient_error <- max(abs(fit$beta - reference))
    prediction_error <- max(abs(features$X %*% (fit$beta - reference)))
    gap_error <- max(abs(fit$primal_dual_gap - fixture$path$diagnostics$primal_dual_gap[seq_len(endpoint)]))
    stopifnot(coefficient_error < 1e-10, prediction_error < 1e-10, gap_error < 1e-10,
      fit$full_grid_size == 40L, fit$solved_grid_size == endpoint, max(abs(fit$primal_dual_gap)) < 1e-7)
    records[[length(records) + 1L]] <- data.table(dataset, endpoint, clusters = features$n,
      variables = ncol(features$X), coefficient_error, prediction_error, gap_error, elapsed_sec = elapsed)
    if (endpoint == 5L) {
      selected <- ep_path(features, "POOLED", 0.75, config, 2L, setting = list(penalty_index = endpoint))
      stopifnot(ncol(selected$beta) == 1L, selected$grid$penalty_index == endpoint,
        max(abs(selected$beta[, 1L] - reference[, endpoint])) < 1e-10)
    }
  }
}
output <- "results/preflight/empirical_selection_20261101/pooled_prefix"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
fwrite(rbindlist(records), file.path(output, "equivalence.csv"))
write_json(list(status = "passed", reference = "saved_complete_native_paths", grid_size = 40L,
  solver = "HiGHS", formal_validation_used = FALSE), file.path(output, "status.json"), auto_unbox = TRUE, pretty = TRUE)
print(rbindlist(records))
