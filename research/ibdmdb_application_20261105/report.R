library(data.table)
library(jsonlite)
library(yaml)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% c(0L, 3L))
config <- read_yaml("config/discovery_ibdmdb_20261105.yml")
implementation <- file.path(config$output_root, "implementation")
for (file in c("prediction_core.R", "ridge_direction.R", "features.R")) source(file.path(implementation, file))
queue <- fread(file.path(config$output_root, "progress.csv"))
stopifnot(nrow(queue) == 35L, all(queue$status == "completed"))
output <- "results/preflight/empirical_selection_20261101/ibdmdb_prediction"
decision_quantile <- config$decision_quantile
if (length(args) == 3L) {
  additional <- fread(file.path(args[1L], "progress.csv"))
  stopifnot(nrow(additional) == 10L, all(additional$status == "completed"))
  queue <- rbindlist(list(queue, additional))
  output <- args[2L]
  decision_quantile <- as.numeric(args[3L])
}
rows <- fread(file.path(config$input_directory, "discovery.csv"))
panel <- ibd_panel(config$input_directory)
records <- checks <- settings <- selected <- list()
for (j in seq_len(nrow(queue))) {
  directory <- queue$directory[j]
  request <- read_json(file.path(directory, "request.json"), simplifyVector = TRUE)
  status <- read_json(file.path(directory, "status.json"), simplifyVector = TRUE)
  fit <- readRDS(file.path(directory, "fit.rds"))
  evaluation <- rows[outer_fold == request$fold]
  stopifnot(status$status == "completed", !request$formal_validation_used,
    identical(fit$recipe$training_strains, sort(unique(rows[outer_fold != request$fold]$strain_id))),
    !length(intersect(fit$recipe$training_strains, evaluation$strain_id)))
  offset <- request$method %in% c("JOINT_OFFSET", "POOLED_OFFSET")
  features <- if (offset) ep_offset_features(evaluation, panel, fit$recipe, ibd_features) else
    ibd_features(evaluation, panel, fit$recipe)
  base <- drop(if (offset) ep_predict_offset(features, fit$beta) else ep_predict(features, fit$beta))
  saved <- fread(file.path(directory, "predictions.csv"))
  error <- 0
  for (k in seq_along(fit$quantiles)) {
    value <- fit$quantiles[k]
    current <- saved[tau == value]
    stopifnot(identical(current$sample_id, evaluation$sample_id), identical(current$strain_id, evaluation$strain_id))
    prediction <- base + fit$calibration[[k]][evaluation$condition]
    error <- max(error, abs(current$y - evaluation$y), abs(current$prediction - prediction),
      abs(current$loss - ep_loss(evaluation$y, prediction, value)))
  }
  stopifnot(is.finite(error), error < 1e-9)
  checks[[j]] <- data.table(fold = request$fold, method = request$method, run = basename(directory), maximum_error = error)
  saved[, `:=`(fold = request$fold, method = request$method)]
  records[[j]] <- saved
  settings[[j]] <- cbind(data.table(fold = request$fold, method = request$method, tau = request$tau,
    runtime_sec = status$runtime_sec), as.data.table(as.data.frame(fit$grid)[fit$choice, , drop = FALSE]))
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
fwrite(crossings, file.path(output, "quantile_crossings.csv"))
write_json(list(status = "completed", phase = "discovery", runs = nrow(queue),
  participants = uniqueN(rows$strain_id), formal_validation_used = FALSE, manuscript_adoption = FALSE,
  uncertainty = "descriptive_participant_bootstrap_of_cross_validated_losses",
  decision_quantile = decision_quantile,
  nonlinear_comparison_required = all(comparison[tau == decision_quantile]$relative_gain > 0)),
  file.path(output, "status.json"), auto_unbox = TRUE, pretty = TRUE)
print(comparison)
