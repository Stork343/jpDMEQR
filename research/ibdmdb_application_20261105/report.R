library(data.table)
library(jsonlite)
library(yaml)
library(grf)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% c(0L, 3L, 4L))
config <- read_yaml("config/discovery_ibdmdb_20261105.yml")
implementation <- file.path(config$output_root, "implementation")
for (file in c("prediction_core.R", "ridge_direction.R", "features.R")) source(file.path(implementation, file))
queue <- fread(file.path(config$output_root, "progress.csv"))
stopifnot(nrow(queue) == 35L, all(queue$status == "completed"))
output <- "results/preflight/empirical_selection_20261101/ibdmdb_prediction"
decision_quantile <- config$decision_quantile
if (length(args) >= 3L) {
  additional <- fread(file.path(args[1L], "progress.csv"))
  stopifnot(nrow(additional) == 10L, all(additional$status == "completed"))
  queue <- rbindlist(list(queue, additional))
  output <- args[2L]
  decision_quantile <- as.numeric(args[3L])
}
if (length(args) == 4L) {
  forest_queue <- fread(file.path(args[4L], "progress.csv"))
  stopifnot(nrow(forest_queue) == 10L, all(forest_queue$status == "completed"))
  queue <- rbindlist(list(queue, forest_queue))
  source(file.path(args[4L], "implementation", "forest.R"))
}
rows <- fread(file.path(config$input_directory, "discovery.csv"))
panel <- ibd_panel(config$input_directory)
records <- checks <- settings <- selected <- forest_checks <- list()
for (j in seq_len(nrow(queue))) {
  directory <- queue$directory[j]
  request <- read_json(file.path(directory, "request.json"), simplifyVector = TRUE)
  status <- read_json(file.path(directory, "status.json"), simplifyVector = TRUE)
  fit <- readRDS(file.path(directory, "fit.rds"))
  evaluation <- rows[outer_fold == request$fold]
  stopifnot(status$status == "completed", !request$formal_validation_used,
    identical(fit$recipe$training_strains, sort(unique(rows[outer_fold != request$fold]$strain_id))),
    !length(intersect(fit$recipe$training_strains, evaluation$strain_id)))
  forest <- request$method %in% c("QRF_RAW", "QRF_OFFSET")
  if (forest) {
    features <- if (fit$offset_model) ep_offset_features(evaluation, panel, fit$recipe, ibd_features) else
      ibd_features(evaluation, panel, fit$recipe$direction_recipe$base_recipe)
    features$sample_id <- evaluation$sample_id
    training <- copy(rows[outer_fold != request$fold])
    roles <- fread(file.path(config$output_root, sprintf("fold%02d", request$fold), "training_roles.csv"))
    stopifnot(identical(roles$sample_id, training$sample_id))
    training[, inner_fold := roles$inner_fold]
    fitting_features <- if (fit$offset_model) ep_offset_features(training, panel, fit$recipe, ibd_features) else
      ibd_features(training, panel, fit$recipe$direction_recipe$base_recipe)
    fitting_features$sample_id <- training$sample_id
    calibration_violation <- 0
    base <- matrix(NA_real_, nrow(evaluation), length(fit$quantiles))
    for (g in unique(fit$settings$forest_index[fit$selected])) {
      model <- fit$models[[as.character(g)]]
      stopifnot(identical(model$sample_ids, training$sample_id),
        model$grf_version == as.character(packageVersion("grf")), identical(model$offset_model, fit$offset_model))
      native <- rs_forest_predict(model, features)
      indices <- which(fit$settings$forest_index[fit$selected] == g)
      base[, indices] <- native[, indices, drop = FALSE]
      modes <- fit$settings$calibration[fit$selected[indices]]
      calibration_predictions <- list()
      for (mode in setdiff(unique(modes), "none")) {
        calibration_predictions[[mode]] <- rs_forest_predict(model, fitting_features, oob = mode == "oob")
      }
      for (k in indices) {
        mode <- fit$settings$calibration[fit$selected[k]]
        calibration <- fit$calibration[[k]]
        if (mode == "none") {
          calibration_violation <- max(calibration_violation, abs(calibration))
        } else {
          residual <- training$y - calibration_predictions[[mode]][, k]
          for (condition in unique(training$condition)) {
            keep <- training$condition == condition
            weight <- fitting_features$weight[keep]
            below <- sum(weight * (residual[keep] < calibration[condition] - 1e-10)) / sum(weight)
            at_or_below <- sum(weight * (residual[keep] <= calibration[condition] + 1e-10)) / sum(weight)
            calibration_violation <- max(calibration_violation, below - fit$quantiles[k], fit$quantiles[k] - at_or_below)
          }
        }
      }
    }
    cv_loss <- matrix(0, nrow(fit$settings), length(fit$quantiles))
    for (fold in sort(unique(training$inner_fold))) {
      validation <- training[inner_fold == fold]
      fitting <- training[inner_fold != fold]
      weight <- 1 / as.numeric(table(validation$strain_id)[validation$strain_id])
      for (g in seq_len(nrow(fit$grid))) {
        inner <- readRDS(file.path(directory, sprintf("inner_%d_forest_%d.rds", fold, g)))
        stopifnot(identical(inner$model$sample_ids, fitting$sample_id),
          !length(intersect(inner$recipe$training_strains, validation$strain_id)))
        for (s in which(fit$settings$forest_index == g)) {
          for (k in seq_along(fit$quantiles)) {
            calibration <- inner$calibration[[paste(s, k, sep = "_")]]
            prediction <- inner$prediction[, k] + calibration[validation$condition]
            cv_loss[s, k] <- cv_loss[s, k] + sum(ep_loss(validation$y, prediction, fit$quantiles[k]) * weight)
          }
        }
        rm(inner)
      }
    }
    cv_error <- max(abs(cv_loss - fit$loss))
    stopifnot(cv_error < 1e-9, calibration_violation < 1e-8,
      identical(as.integer(apply(cv_loss, 2L, which.min)), as.integer(fit$selected)))
    forest_checks[[length(forest_checks) + 1L]] <- data.table(fold = request$fold,
      method = request$method, maximum_cv_loss_error = cv_error, maximum_calibration_violation = calibration_violation)
    settings[[j]] <- rbindlist(lapply(seq_along(fit$quantiles), function(k) {
      current <- fit$settings[fit$selected[k], ]
      cbind(data.table(fold = request$fold, method = request$method, tau = fit$quantiles[k],
        runtime_sec = status$runtime_sec), current, as.data.table(fit$grid[current$forest_index, , drop = FALSE]))
    }))
  } else {
    offset <- request$method %in% c("JOINT_OFFSET", "POOLED_OFFSET")
    features <- if (offset) ep_offset_features(evaluation, panel, fit$recipe, ibd_features) else
      ibd_features(evaluation, panel, fit$recipe)
    base <- drop(if (offset) ep_predict_offset(features, fit$beta) else ep_predict(features, fit$beta))
    settings[[j]] <- cbind(data.table(fold = request$fold, method = request$method, tau = request$tau,
      runtime_sec = status$runtime_sec), as.data.table(as.data.frame(fit$grid)[fit$choice, , drop = FALSE]))
  }
  saved <- fread(file.path(directory, "predictions.csv"))
  error <- 0
  for (k in seq_along(fit$quantiles)) {
    value <- fit$quantiles[k]
    current <- saved[tau == value]
    stopifnot(identical(current$sample_id, evaluation$sample_id), identical(current$strain_id, evaluation$strain_id))
    prediction <- (if (forest) base[, k] else base) + fit$calibration[[k]][evaluation$condition]
    error <- max(error, abs(current$y - evaluation$y), abs(current$prediction - prediction),
      abs(current$loss - ep_loss(evaluation$y, prediction, value)))
  }
  stopifnot(is.finite(error), error < 1e-9)
  checks[[j]] <- data.table(fold = request$fold, method = request$method, run = basename(directory), maximum_error = error)
  saved[, `:=`(fold = request$fold, method = request$method)]
  records[[j]] <- saved
  if (request$method == "JOINT_OFFSET") {
    active <- which(abs(fit$beta) > 1e-10)
    selected[[length(selected) + 1L]] <- data.table(fold = request$fold, tau = request$tau,
      variable = names(fit$beta)[active], coefficient_scaled = fit$beta[active])
  }
}
predictions <- rbindlist(records)
stopifnot(!anyDuplicated(predictions[, .(sample_id, method, tau)]))
losses <- predictions[, .(loss = mean(loss)), by = .(strain_id, fold, method, tau)]
counts <- losses[, .N, by = .(method, tau)]
stopifnot(all(counts$N == uniqueN(rows$strain_id)))
summary <- losses[, .(participants = .N, loss = mean(loss), se = sd(loss) / sqrt(.N)), by = .(method, tau)]
set.seed(20261105L)
weights <- rmultinom(10000L, uniqueN(rows$strain_id), rep(1 / uniqueN(rows$strain_id), uniqueN(rows$strain_id))) / uniqueN(rows$strain_id)
comparisons <- list()
for (value in sort(unique(losses$tau))) {
  wide <- dcast(losses[tau == value], strain_id ~ method, value.var = "loss")
  stopifnot(!anyNA(wide), identical(wide$strain_id, sort(unique(rows$strain_id))))
  for (method_name in setdiff(names(wide), c("strain_id", "JOINT_OFFSET"))) {
    reference <- wide[[method_name]]
    difference <- reference - wide$JOINT_OFFSET
    bootstrap <- drop(crossprod(weights, difference)) / drop(crossprod(weights, reference))
    comparisons[[length(comparisons) + 1L]] <- data.table(tau = value, comparator = method_name,
      participants = nrow(wide), joint_loss = mean(wide$JOINT_OFFSET), comparator_loss = mean(reference),
      relative_gain = mean(difference) / mean(reference), lower = quantile(bootstrap, 0.025),
      upper = quantile(bootstrap, 0.975), phase = "discovery", formal_validation_used = FALSE, manuscript_adoption = FALSE)
  }
}
comparison <- rbindlist(comparisons)
crossings <- predictions[, {
  available <- sort(unique(tau))
  if (identical(available, c(0.25, 0.5, 0.75))) {
    wide <- dcast(.SD, sample_id ~ tau, value.var = "prediction")
    gaps <- cbind(wide[["0.25"]] - wide[["0.5"]], wide[["0.5"]] - wide[["0.75"]])
    list(quantiles = 3L, samples = nrow(wide), crossed_samples = sum(apply(gaps, 1L, max) > 1e-12), maximum_crossing = max(0, gaps))
  } else list(quantiles = length(available), samples = uniqueN(sample_id), crossed_samples = NA_integer_, maximum_crossing = NA_real_)
}, by = method]
dir.create(output, recursive = TRUE, showWarnings = FALSE)
fwrite(summary, file.path(output, "summary.csv"))
fwrite(comparison, file.path(output, "comparisons.csv"))
fwrite(rbindlist(checks), file.path(output, "artifact_checks.csv"))
fwrite(rbindlist(settings, fill = TRUE), file.path(output, "settings.csv"))
fwrite(rbindlist(selected), file.path(output, "selected_coefficients.csv"))
if (length(forest_checks)) fwrite(rbindlist(forest_checks), file.path(output, "forest_cv_checks.csv"))
fwrite(crossings, file.path(output, "quantile_crossings.csv"))
write_json(list(status = "completed", phase = "discovery", runs = nrow(queue),
  participants = uniqueN(rows$strain_id), formal_validation_used = FALSE, manuscript_adoption = FALSE,
  uncertainty = "descriptive_participant_bootstrap_of_cross_validated_losses",
  decision_quantile = decision_quantile,
  nonlinear_comparison_required = all(comparison[tau == decision_quantile]$relative_gain > 0)),
  file.path(output, "status.json"), auto_unbox = TRUE, pretty = TRUE)
print(comparison)
