library(yaml)
library(jsonlite)
library(reticulate)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
config <- read_yaml(args[1L])
study <- read_yaml(args[2L])
fold <- as.integer(args[3L])
implementation <- file.path(config$output_root, "implementation")
source(file.path(implementation, "engine.R"))
source(file.path(implementation, "panel.R"))
source(file.path(implementation, "prediction_core.R"))
use_python(file.path(getwd(), "results/scratch/cq_solver_env/bin/python"), required = TRUE)
Sys.setenv(NUMBA_CACHE_DIR = normalizePath("results/scratch/runtime", mustWork = TRUE))
source_python(file.path(implementation, "smooth_path.py"))
stopifnot(config$smoothing_exponent == 0.2, config$numerical_ridge == "inverse_training_patients")
outer <- file.path(study$output_root, paste0("outer_", fold))
rows <- fread(file.path(outer, "training.csv"))
evaluation <- fread(file.path(outer, "evaluation.csv"))
recipes <- readRDS(file.path(outer, "recipes.rds"))
stopifnot(all(rows$role == "training"), all(evaluation$role == "discovery_evaluation"),
  !length(intersect(rows$strain_id, evaluation$strain_id)))
panel <- cp_panel(study$input_directory)
root <- file.path(config$output_root, study$study_id, paste0("outer_", fold))
dir.create(root, recursive = TRUE, showWarnings = FALSE)

smooth_path <- function(features, tau, alpha, indices = seq_len(config$path_length)) {
  lambda <- cq_lambda(features, tau, alpha, config)[indices]
  clinical <- which(features$penalty_factor == 0)
  native <- rq.fit.br(features$X[, clinical, drop = FALSE] * features$weight,
    features$y * features$weight, tau = tau)
  initial <- numeric(ncol(features$X))
  initial[clinical] <- native$coefficients
  fitted <- fit_smooth_path(features$X, features$y, features$weight, features$penalty_factor,
    lambda, alpha, tau, features$n, initial, config$solver_tolerance)
  stopifnot(max(fitted$diagnostics[, 4L]) < config$maximum_kkt)
  rownames(fitted$beta) <- colnames(features$X)
  c(fitted, list(lambda = lambda, alpha = alpha))
}

for (tau in unlist(config$quantiles)) {
  directory <- file.path(root, sprintf("tau%03d", as.integer(tau * 100)))
  if (file.exists(file.path(directory, "status.json"))) {
    stopifnot(read_json(file.path(directory, "status.json"))$status == "completed")
    next
  }
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  started <- proc.time()[["elapsed"]]
  grid <- CJ(alpha = unlist(config$alpha), index = seq_len(config$path_length), sorted = FALSE)
  loss <- numeric(nrow(grid))
  oof <- matrix(NA_real_, nrow(rows), nrow(grid))
  for (inner in sort(unique(rows$inner_fold))) {
    recipe <- recipes[[paste0("inner_", inner)]]$direction_recipe$base_recipe
    training_rows <- rows[inner_fold != inner]
    testing_rows <- rows[inner_fold == inner]
    stopifnot(identical(sort(unique(training_rows$strain_id)), recipe$training_strains))
    training <- cq_features(training_rows, panel, recipe)
    testing <- cq_features(testing_rows, panel, recipe, colnames(training$X))
    for (alpha in unlist(config$alpha)) {
      path_file <- file.path(directory, sprintf("inner_%d_alpha_%s.rds", inner, alpha))
      if (file.exists(path_file)) {
        saved <- readRDS(path_file)
        stopifnot(identical(saved$columns, colnames(training$X)), saved$tau == tau, saved$path$alpha == alpha)
        path <- saved$path
      } else {
        cat(study$study_id, fold, tau, "inner", inner, "alpha", alpha, "\n")
        flush.console()
        path <- smooth_path(training, tau, alpha)
        saveRDS(list(path = path, columns = colnames(training$X), recipe = recipe, tau = tau), paste0(path_file, ".partial"))
        stopifnot(file.rename(paste0(path_file, ".partial"), path_file))
      }
      predictions <- ep_predict(testing, path$beta)
      indices <- which(grid$alpha == alpha)
      loss[indices] <- loss[indices] + colSums(ep_loss(testing_rows$y, predictions, tau) * testing$weight)
      oof[rows$inner_fold == inner, indices] <- predictions
    }
  }
  stopifnot(all(is.finite(oof)))
  choice <- which.min(loss)
  setting <- grid[choice]
  recipe <- recipes$training$direction_recipe$base_recipe
  training <- cq_features(rows, panel, recipe)
  testing <- cq_features(evaluation, panel, recipe, colnames(training$X))
  final <- smooth_path(training, tau, setting$alpha, seq_len(setting$index))
  beta <- final$beta[, setting$index]
  prediction <- drop(ep_predict(testing, beta))
  fitted <- list(beta = beta, recipe = recipe, columns = colnames(training$X), grid = grid,
    inner_loss = loss / uniqueN(rows$strain_id), choice = choice, oof = oof[, choice],
    alpha = setting$alpha, lambda = final$lambda[setting$index], tau = tau, prediction = prediction,
    bandwidth = final$bandwidth, diagnostics = final$diagnostics[setting$index, ],
    clinical_columns = names(beta)[training$penalty_factor == 0])
  saveRDS(fitted, file.path(directory, "fit.rds"))
  records <- data.table(study_id = study$study_id, fold, method = config$prediction_method, tau,
    sample_id = evaluation$sample_id, strain_id = evaluation$strain_id, y = evaluation$y, prediction)
  fwrite(records, file.path(directory, "prediction.csv"))
  diagnostic <- data.table(study_id = study$study_id, fold, tau, alpha = setting$alpha,
    lambda = fitted$lambda, selected_index = setting$index, bandwidth = final$bandwidth,
    nonzero_proteins = sum(abs(beta) > 1e-7 & startsWith(names(beta), "P_")),
    clinical_columns = length(fitted$clinical_columns), maximum_kkt = fitted$diagnostics[4L],
    objective = fitted$diagnostics[5L], runtime_sec = proc.time()[["elapsed"]] - started)
  fwrite(diagnostic, file.path(directory, "diagnostics.csv"))
  write_json(list(status = "completed", formal_validation_used = FALSE, manuscript_adoption = FALSE),
    file.path(directory, "status.json"), auto_unbox = TRUE, pretty = TRUE)
  cat(study$study_id, fold, tau, "completed", diagnostic$runtime_sec, "seconds\n")
  flush.console()
  gc()
}
writeLines(trimws(capture.output(sessionInfo()), which = "right"), file.path(root, "sessionInfo.txt"))
