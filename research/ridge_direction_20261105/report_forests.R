library(data.table)
library(jsonlite)
library(yaml)
library(grf)

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
assay_name <- args[1L]
config <- read_yaml("config/development_ridge_offset_forest_20261105.yml")
stopifnot(assay_name %in% config$assays)
implementation <- file.path(config$output_root, "implementation")
source(file.path(implementation, "prediction_core.R"))
source(file.path(implementation, "ridge_direction.R"))
source(file.path(implementation, "contrasts.R"))
source(file.path(implementation, "features.R"))
source(file.path(implementation, "forest.R"))

queue <- fread(file.path(config$output_root, "progress.csv"))[trait == assay_name]
stopifnot(nrow(queue) == 2L, all(queue$status == "completed"))
linear <- fread(file.path(config$source_recipes, "progress.csv"))[trait == assay_name]
prior <- data.table(trait = assay_name, method = "ENET", status = "completed",
  directory = file.path("results/raw/application/yeast_proteome_20261105/discovery", assay_name, "ENET_weighted_cv"))
queue <- rbindlist(list(queue, linear, prior), fill = TRUE)
stopifnot(all(queue$status == "completed"), !anyDuplicated(queue$directory))
rows <- fread(file.path(config$input_directory, assay_name, "discovery.csv"))
stopifnot(all(rows$role == "discovery"))
training <- rows[exploration_role == "training"]
evaluation <- rows[exploration_role == "evaluation"]
stopifnot(!length(intersect(training$strain_id, evaluation$strain_id)), all(evaluation$y >= 0))
panel <- proteome_panel(config$input_directory)
recipes <- readRDS(file.path(config$source_recipes, assay_name, "recipes.rds"))
records <- checks <- settings <- cv_checks <- list()

for (j in seq_len(nrow(queue))) {
  directory <- queue$directory[j]
  request <- read_json(file.path(directory, "request.json"), simplifyVector = TRUE)
  status <- read_json(file.path(directory, "status.json"), simplifyVector = TRUE)
  stopifnot(status$status == "completed", !request$formal_validation_used)
  fit <- readRDS(file.path(directory, "fit.rds"))
  stopifnot(identical(fit$recipe$training_strains, sort(unique(training$strain_id))))
  saved <- fread(file.path(directory, "predictions.csv"))
  method_name <- request$method
  forest <- method_name %in% c("QRF_RAW", "QRF_OFFSET")
  maximum_error <- 0
  if (forest) {
    feature_function <- function(current, recipe) {
      result <- if (fit$offset_model) ep_offset_features(current, panel, recipe, proteome_features) else
        proteome_features(current, panel, recipe$direction_recipe$base_recipe)
      result$sample_id <- current$sample_id
      result
    }
    features <- feature_function(evaluation, fit$recipe)
    base <- matrix(NA_real_, nrow(evaluation), length(fit$quantiles))
    for (g in unique(fit$settings$forest_index[fit$selected])) {
      model <- fit$models[[as.character(g)]]
      stopifnot(identical(model$sample_ids, training$sample_id),
                model$grf_version == as.character(packageVersion("grf")),
                identical(model$offset_model, fit$offset_model))
      native_prediction <- rs_forest_predict(model, features)
      selected_quantiles <- which(fit$settings$forest_index[fit$selected] == g)
      base[, selected_quantiles] <- native_prediction[, selected_quantiles, drop = FALSE]
    }
    reconstructed_loss <- matrix(0, nrow(fit$settings), length(fit$quantiles))
    for (fold in sort(unique(training$inner_fold))) {
      validation <- training[inner_fold == fold]
      fitting <- training[inner_fold != fold]
      weight <- 1 / as.numeric(table(validation$strain_id)[validation$strain_id])
      for (g in seq_len(nrow(fit$grid))) {
        inner <- readRDS(file.path(directory, sprintf("inner_%d_forest_%d.rds", fold, g)))
        stopifnot(identical(inner$model$sample_ids, fitting$sample_id),
                  !length(intersect(inner$recipe$training_strains, validation$strain_id)),
                  identical(dim(inner$prediction), c(nrow(validation), length(fit$quantiles))))
        for (s in which(fit$settings$forest_index == g)) {
          for (k in seq_along(fit$quantiles)) {
            calibration <- inner$calibration[[paste(s, k, sep = "_")]]
            prediction <- inner$prediction[, k] + calibration[validation$condition]
            reconstructed_loss[s, k] <- reconstructed_loss[s, k] +
              sum(ep_loss(validation$y, prediction, fit$quantiles[k]) * weight)
          }
        }
        rm(inner)
      }
    }
    cv_error <- max(abs(reconstructed_loss - fit$loss))
    stopifnot(cv_error < 1e-9, identical(as.integer(apply(reconstructed_loss, 2L, which.min)), as.integer(fit$selected)))
    cv_checks[[length(cv_checks) + 1L]] <- data.table(assay = assay_name, method = method_name,
      inner_folds = uniqueN(training$inner_fold), forest_settings = nrow(fit$grid),
      calibration_settings = nrow(fit$settings), maximum_error = cv_error)
    for (k in seq_along(fit$quantiles)) {
      selected <- fit$settings[fit$selected[k], ]
      settings[[length(settings) + 1L]] <- cbind(data.table(assay = assay_name, method = method_name,
        tau = fit$quantiles[k], runtime_sec = status$runtime_sec), selected,
        as.data.table(fit$grid[selected$forest_index, , drop = FALSE]))
    }
  } else {
    features <- if (method_name %in% c("JOINT_OFFSET", "POOLED_OFFSET"))
      ep_offset_features(evaluation, panel, fit$recipe, proteome_features) else
      proteome_features(evaluation, panel, fit$recipe)
    prediction <- if (method_name %in% c("JOINT_OFFSET", "POOLED_OFFSET"))
      ep_predict_offset(features, fit$beta) else ep_predict(features, fit$beta)
    base <- matrix(prediction, nrow(evaluation), length(fit$quantiles))
  }
  for (k in seq_along(fit$quantiles)) {
    value <- fit$quantiles[k]
    current <- saved[tau == value]
    stopifnot(identical(current$sample_id, evaluation$sample_id),
              identical(current$strain_id, evaluation$strain_id), identical(current$condition, evaluation$condition))
    prediction <- base[, k] + fit$calibration[[k]][evaluation$condition]
    maximum_error <- max(maximum_error, abs(current$y - evaluation$y),
      abs(current$prediction - prediction), abs(current$loss - ep_loss(evaluation$y, prediction, value)))
  }
  stopifnot(is.finite(maximum_error), maximum_error < 1e-9,
            nrow(saved) == nrow(evaluation) * length(fit$quantiles))
  checks[[length(checks) + 1L]] <- data.table(assay = assay_name, method = method_name,
    run = basename(directory), rows = nrow(saved), maximum_error = maximum_error)
  saved[, `:=`(assay = assay_name, method = method_name)]
  records[[j]] <- saved
  rm(features, fit)
  gc()
}

predictions <- rbindlist(records)
stopifnot(!anyDuplicated(predictions[, .(assay, method, sample_id, tau)]))
projected <- copy(predictions)
projected[, prediction := pmax(prediction, 0)]
projected[, loss := ep_loss(y, prediction, tau)]
stopifnot(all(projected$loss <= predictions$loss + 1e-12))
predictions[, scoring := "registered"]
projected[, scoring := "nonnegative_sensitivity"]
values <- rbindlist(list(predictions, projected))
losses <- values[, .(loss = mean(loss)), by = .(assay, scoring, method, strain_id, tau)]
summary <- losses[, .(strains = .N, loss = mean(loss), se = sd(loss) / sqrt(.N)),
  by = .(assay, scoring, method, tau)]
set.seed(20261105L)
bootstrap_weights <- rmultinom(10000L, uniqueN(evaluation$strain_id),
  rep(1 / uniqueN(evaluation$strain_id), uniqueN(evaluation$strain_id))) / uniqueN(evaluation$strain_id)
comparisons <- list()
for (score_name in unique(losses$scoring)) {
  for (value in sort(unique(losses$tau))) {
    wide <- dcast(losses[scoring == score_name & tau == value], strain_id ~ method, value.var = "loss")
    stopifnot(!anyNA(wide), identical(wide$strain_id, sort(unique(evaluation$strain_id))))
    for (method_name in setdiff(names(wide), c("strain_id", "JOINT_OFFSET"))) {
      reference <- wide[[method_name]]
      difference <- reference - wide$JOINT_OFFSET
      bootstrap <- drop(crossprod(bootstrap_weights, difference)) / drop(crossprod(bootstrap_weights, reference))
      comparisons[[length(comparisons) + 1L]] <- data.table(assay = assay_name, scoring = score_name,
        tau = value, target_method = "JOINT_OFFSET", comparator = method_name, strains = nrow(wide),
        joint_loss = mean(wide$JOINT_OFFSET), comparator_loss = mean(reference),
        relative_gain = mean(difference) / mean(reference), lower = quantile(bootstrap, 0.025),
        upper = quantile(bootstrap, 0.975), phase = "method_development", manuscript_adoption = FALSE)
    }
  }
}
comparison <- rbindlist(comparisons)
crossings <- predictions[, {
  available <- sort(unique(tau))
  if (identical(available, c(0.25, 0.5, 0.75))) {
    wide <- dcast(.SD, sample_id ~ tau, value.var = "prediction")
    gaps <- cbind(wide[["0.25"]] - wide[["0.5"]], wide[["0.5"]] - wide[["0.75"]])
    list(quantiles = length(available), records = nrow(wide),
         crossed_records = sum(apply(gaps, 1L, max) > 1e-12), maximum_crossing = max(0, gaps))
  } else list(quantiles = length(available), records = uniqueN(sample_id), crossed_records = NA_integer_, maximum_crossing = NA_real_)
}, by = .(assay, method)]
condition_metrics <- values[, .(records = .N, strains = uniqueN(strain_id), loss = mean(loss),
  fraction_below = mean(y < prediction), fraction_at_or_below = mean(y <= prediction),
  negative_predictions = sum(prediction < 0)), by = .(assay, scoring, method, tau, condition)]
decision <- comparison[scoring == "registered" & tau == config$decision_quantile,
  .(comparators = .N, smallest_relative_gain = min(relative_gain), smallest_interval_lower = min(lower),
    point_gain_all_comparators = all(relative_gain > 0), positive_interval_all_comparators = all(lower > 0)), by = assay]
stopifnot(decision$comparators == 5L)
decision[, `:=`(formal_validation_used = FALSE, manuscript_adoption = FALSE)]
output <- file.path("results/preflight/empirical_selection_20261101/ridge_offset_forests", assay_name)
dir.create(output, recursive = TRUE, showWarnings = FALSE)
fwrite(summary, file.path(output, "summary.csv"))
fwrite(comparison, file.path(output, "comparisons.csv"))
fwrite(rbindlist(checks), file.path(output, "artifact_checks.csv"))
fwrite(rbindlist(cv_checks), file.path(output, "cv_checks.csv"))
fwrite(rbindlist(settings), file.path(output, "forest_settings.csv"))
fwrite(crossings, file.path(output, "quantile_crossings.csv"))
fwrite(condition_metrics, file.path(output, "condition_metrics.csv"))
fwrite(decision, file.path(output, "screening_decision.csv"))
write_json(list(status = "completed", phase = "method_development", assay = assay_name,
  completed_runs = nrow(queue), formal_validation_used = FALSE, manuscript_adoption = FALSE,
  interval_type = "exploratory_paired_strain_bootstrap_95_percent", bootstrap_repetitions = 10000L,
  primary_scoring = "registered", sensitivity = "nonnegative_projection_of_all_predictions"),
  file.path(output, "status.json"), auto_unbox = TRUE, pretty = TRUE)
print(comparison[scoring == "registered" & tau == config$decision_quantile])
print(decision)
