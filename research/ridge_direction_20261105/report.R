library(data.table)
library(jsonlite)
library(yaml)
source("research/yeast_stress_application_20261104/features.R")
source("research/yeast_proteome_application_20261105/features.R")
source("research/empirical_selection_20261101/prediction_core.R")
source("research/empirical_selection_20261101/ridge_direction.R")

config <- read_yaml("config/development_ridge_direction_20261105.yml")
panel <- proteome_panel(config$input_directory)
evaluation <- fread(file.path(config$input_directory, "NaCl", "discovery.csv"))[exploration_role == "evaluation"]
training <- fread(file.path(config$input_directory, "NaCl", "discovery.csv"))[exploration_role == "training"]
first <- fread("results/raw/application/ridge_direction_20261105/progress.csv")
second <- fread("results/raw/application/ridge_offset_20261105/progress.csv")
stopifnot(all(first$status == "completed"), all(second$status == "completed"))
directories <- c(first$directory, second$directory,
  "results/raw/application/ridge_direction_20261105/NaCl/ADJUSTMENT_MEAN_cv",
  "results/raw/application/yeast_proteome_20261105/discovery/NaCl/JOINT_prediction_tau025",
  "results/raw/application/yeast_proteome_20261105/discovery/NaCl/ENET_weighted_cv",
  file.path("results/raw/application/calibrated_profile_20261105",
    c("JOINT_CALIBRATED", "DIRECTION_CALIBRATED", "OFFSET_CALIBRATED")))
records <- checks <- settings <- list()
for (j in seq_along(directories)) {
  directory <- directories[j]
  status <- read_json(file.path(directory, "status.json"), simplifyVector = TRUE)
  request <- read_json(file.path(directory, "request.json"), simplifyVector = TRUE)
  stopifnot(status$status == "completed", !request$formal_validation_used)
  method <- request$method
  fit <- readRDS(file.path(directory, "fit.rds"))
  if (endsWith(method, "OFFSET") || method == "OFFSET_CALIBRATED") {
    features <- ep_offset_features(evaluation, panel, fit$recipe, proteome_features)
    predicted <- drop(ep_predict_offset(features, fit$beta))
  } else if (endsWith(method, "DIRECTION") || method == "DIRECTION_CALIBRATED") {
    features <- ep_direction_features(evaluation, panel, fit$recipe, proteome_features)
    predicted <- drop(ep_predict(features, fit$beta))
  } else {
    features <- proteome_features(evaluation, panel, fit$recipe)
    predicted <- drop(ep_predict(features, fit$beta))
  }
  saved <- fread(file.path(directory, "predictions.csv"))
  error <- 0
  for (k in seq_along(fit$quantiles)) {
    value <- fit$quantiles[k]
    current <- saved[tau == value]
    stopifnot(identical(current$sample_id, evaluation$sample_id))
    reconstructed <- predicted + fit$calibration[[k]][evaluation$condition]
    error <- max(error, abs(current$y - evaluation$y), abs(current$prediction - reconstructed),
                 abs(current$loss - ep_loss(evaluation$y, reconstructed, value)))
  }
  stopifnot(error < 1e-9)
  calibration_error <- NA_real_
  if (endsWith(method, "CALIBRATED")) {
    if (method == "OFFSET_CALIBRATED") {
      fitting <- ep_offset_features(training, panel, fit$recipe, proteome_features)
      fitted_prediction <- drop(ep_predict_offset(fitting, fit$beta))
    } else if (method == "DIRECTION_CALIBRATED") {
      fitting <- ep_direction_features(training, panel, fit$recipe, proteome_features)
      fitted_prediction <- drop(ep_predict(fitting, fit$beta))
    } else {
      fitting <- proteome_features(training, panel, fit$recipe)
      fitted_prediction <- drop(ep_predict(fitting, fit$beta))
    }
    residual <- training$y - fitted_prediction
    weight <- 1 / as.numeric(table(training$strain_id)[training$strain_id])
    calibration_error <- 0
    for (condition_name in unique(training$condition)) {
      indices <- which(training$condition == condition_name)
      alpha <- fit$calibration[[1L]][condition_name]
      target <- fit$quantiles[1L] * sum(weight[indices])
      lower <- sum(weight[indices][residual[indices] < alpha - 1e-10])
      upper <- sum(weight[indices][residual[indices] <= alpha + 1e-10])
      calibration_error <- max(calibration_error, lower - target, target - upper)
    }
    stopifnot(calibration_error < 1e-8)
  }
  checks[[j]] <- data.table(method = method, run = basename(directory), maximum_error = error,
                            training_calibration_kkt_error = calibration_error)
  loss <- saved[, .(loss = mean(loss)), by = .(strain_id, tau)]
  loss[, method := method]
  records[[j]] <- loss
  if (!is.null(fit$grid)) {
    requested_tau <- if (length(fit$quantiles) == 1L) fit$quantiles else NA_real_
    settings[[j]] <- cbind(data.table(method = method, tau = requested_tau, runtime_sec = status$runtime_sec),
      as.data.table(as.data.frame(fit$grid)[fit$choice, , drop = FALSE]))
  }
}
losses <- rbindlist(records)
stopifnot(!anyDuplicated(losses[, .(strain_id, tau, method)]))
summary <- losses[, .(strains = .N, loss = mean(loss), se = sd(loss) / sqrt(.N)), by = .(method, tau)]
set.seed(20261105L)
comparisons <- list()
for (value in sort(unique(losses$tau))) {
  wide <- dcast(losses[tau == value], strain_id ~ method, value.var = "loss")
  stopifnot(!anyNA(wide), nrow(wide) == 143L)
  weights <- rmultinom(10000L, nrow(wide), rep(1 / nrow(wide), nrow(wide))) / nrow(wide)
  for (model in c("JOINT_DIRECTION", "JOINT_OFFSET", "JOINT_CALIBRATED", "DIRECTION_CALIBRATED", "OFFSET_CALIBRATED")) {
    if (!(model %in% names(wide))) next
    for (reference_method in setdiff(names(wide), c("strain_id", model))) {
      reference <- wide[[reference_method]]
      difference <- reference - wide[[model]]
      bootstrap <- drop(crossprod(weights, difference)) / drop(crossprod(weights, reference))
      comparisons[[length(comparisons) + 1L]] <- data.table(method = model, tau = value,
        comparator = reference_method, loss = mean(wide[[model]]), comparator_loss = mean(reference),
        relative_gain = mean(difference) / mean(reference), lower = quantile(bootstrap, 0.025),
        upper = quantile(bootstrap, 0.975), phase = "method_development", manuscript_adoption = FALSE)
    }
  }
}
output <- "results/preflight/empirical_selection_20261101/ridge_direction_comparisons"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
comparison <- rbindlist(comparisons)
fwrite(summary, file.path(output, "summary.csv"))
fwrite(comparison, file.path(output, "comparisons.csv"))
fwrite(rbindlist(checks), file.path(output, "artifact_checks.csv"))
fwrite(rbindlist(settings, fill = TRUE), file.path(output, "settings.csv"))
write_json(list(stage = "method_development", completed_runs = length(directories),
                formal_validation_used = FALSE, manuscript_adoption = FALSE),
           file.path(output, "status.json"), auto_unbox = TRUE, pretty = TRUE)
print(summary)
print(comparison[method == "JOINT_OFFSET" & tau == 0.25])
