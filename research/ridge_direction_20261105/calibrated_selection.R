library(data.table)
library(jsonlite)
library(yaml)
source("research/yeast_stress_application_20261104/features.R")
source("research/yeast_proteome_application_20261105/features.R")
source("research/empirical_selection_20261101/prediction_core.R")
source("research/empirical_selection_20261101/ridge_direction.R")
source("research/score_validity_20261026/load_joint.R")
source("research/joint_profile_manuscript_20261026/conic_search.R")

base <- read_yaml("config/discovery_yeast_proteome_20261105.yml")
base$conditions <- base$assays$NaCl
panel <- proteome_panel(base$input_directory)
rows <- fread(file.path(base$input_directory, "NaCl", "discovery.csv"))
training <- rows[exploration_role == "training"]
evaluation <- rows[exploration_role == "evaluation"]
tau <- 0.25
sources <- c(
  JOINT_CALIBRATED = "results/raw/application/yeast_proteome_20261105/discovery/NaCl/JOINT_prediction_tau025",
  DIRECTION_CALIBRATED = "results/raw/application/ridge_direction_20261105/NaCl/JOINT_DIRECTION_tau025",
  OFFSET_CALIBRATED = "results/raw/application/ridge_offset_20261105/NaCl/JOINT_OFFSET_tau025")
root <- "results/raw/application/calibrated_profile_20261105"
if (dir.exists(root)) stop("开发实验目录已经存在")
dir.create(root, recursive = TRUE)
summaries <- list()
for (method in names(sources)) {
  source_directory <- sources[[method]]
  directory <- file.path(root, method)
  dir.create(directory)
  feature_function <- if (method == "JOINT_CALIBRATED") {
    proteome_features
  } else if (method == "DIRECTION_CALIBRATED") {
    function(rows, panel, recipe) ep_direction_features(rows, panel, recipe, proteome_features)
  } else {
    function(rows, panel, recipe) ep_offset_features(rows, panel, recipe, proteome_features)
  }
  prediction_function <- if (method == "OFFSET_CALIBRATED") ep_predict_offset else ep_predict
  started <- proc.time()[["elapsed"]]
  loss <- NULL
  for (fold in sort(unique(training$inner_fold))) {
    saved <- readRDS(file.path(source_directory, paste0("inner_", fold, ".rds")))
    fitting_rows <- training[inner_fold != fold]
    validation_rows <- training[inner_fold == fold]
    fitting <- feature_function(fitting_rows, panel, saved$recipe)
    predicted_training <- as.matrix(prediction_function(fitting, saved$path$beta))
    predicted_validation <- saved$prediction
    stopifnot(nrow(predicted_validation) == nrow(validation_rows))
    if (is.null(loss)) {
      loss <- numeric(ncol(predicted_validation))
      grid <- saved$path$grid
    }
    stopifnot(identical(as.list(grid), as.list(saved$path$grid)))
    validation_weight <- 1 / as.numeric(table(validation_rows$strain_id)[validation_rows$strain_id])
    for (j in seq_len(ncol(predicted_validation))) {
      calibration <- ep_calibrate(fitting_rows, predicted_training[, j], tau, "strain_id", "condition")
      adjusted <- predicted_validation[, j] + calibration[validation_rows$condition]
      loss[j] <- loss[j] + sum(ep_loss(validation_rows$y, adjusted, tau) * validation_weight)
    }
  }
  choice <- which.min(loss)
  original <- readRDS(file.path(source_directory, "fit.rds"))
  recipe <- original$recipe
  fitting <- feature_function(training, panel, recipe)
  final <- ep_path(fitting, "JOINT", tau, base, max(grid$size), as.list(grid[choice]))
  beta <- as.numeric(final$beta[, 1L])
  names(beta) <- colnames(fitting$X)
  calibration <- ep_calibrate(training, drop(prediction_function(fitting, beta)), tau, "strain_id", "condition")
  test <- feature_function(evaluation, panel, recipe)
  predicted <- drop(prediction_function(test, beta)) + calibration[evaluation$condition]
  predictions <- data.table(sample_id = evaluation$sample_id, strain_id = evaluation$strain_id,
    condition = evaluation$condition, y = evaluation$y, tau = tau, prediction = predicted,
    loss = ep_loss(evaluation$y, predicted, tau))
  fit <- list(beta = beta, recipe = recipe, prediction = matrix(predicted, ncol = 1L), quantiles = tau,
              calibration = list(calibration), choice = choice, grid = grid, loss = loss,
              final_diagnostics = final$diagnostics)
  saveRDS(fit, file.path(directory, "fit.rds"))
  fwrite(predictions, file.path(directory, "predictions.csv"))
  losses <- predictions[, .(loss = mean(loss)), by = .(strain_id, tau)]
  fwrite(losses, file.path(directory, "strain_loss.csv"))
  elapsed <- proc.time()[["elapsed"]] - started
  summary <- data.table(method = method, tau = tau, strains = nrow(losses), loss = mean(losses$loss),
    se = sd(losses$loss) / sqrt(nrow(losses)), runtime_sec = elapsed,
    Lambda = grid$Lambda[choice], size = grid$size[choice], phase = "method_development")
  fwrite(summary, file.path(directory, "summary.csv"))
  write_json(list(method = method, source_directory = source_directory,
    calibration = "training_residual_per_condition_in_each_fit_and_final_fit",
    selection_loss = "calibrated_validation_pinball", code_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
    formal_validation_used = FALSE, manuscript_adoption = FALSE),
    file.path(directory, "request.json"), auto_unbox = TRUE, pretty = TRUE)
  write_json(list(status = "completed", runtime_sec = elapsed, formal_validation_used = FALSE,
                  manuscript_adoption = FALSE), file.path(directory, "status.json"), auto_unbox = TRUE, pretty = TRUE)
  summaries[[method]] <- summary
  print(summary)
}
fwrite(rbindlist(summaries), file.path(root, "summary.csv"))
