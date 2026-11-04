library(data.table)
library(jsonlite)
library(yaml)
library(BEDMatrix)
source("research/yeast_1011_application_20261104/features.R")
source("research/empirical_selection_20261101/prediction_core.R")

config <- read_yaml("config/discovery_yeast_1011_20261104.yml")
queue <- fread("results/raw/application/yeast_1011_20261104/screening_20261104/progress.csv")
stopifnot(nrow(queue) == 5L, all(queue$status == "completed"))
panel <- y1011_panel(config$input_directory)
evaluation <- fread(file.path(config$input_directory, "discovery.csv"))[exploration_role == "evaluation"]
stopifnot(all(evaluation$role == "discovery"))
records <- checks <- selections <- condition_records <- list()
for (directory in queue$directory) {
  request <- read_json(file.path(directory, "request.json"), simplifyVector = TRUE)
  status <- read_json(file.path(directory, "status.json"), simplifyVector = TRUE)
  stopifnot(status$status == "completed", !request$formal_validation_used)
  fit <- readRDS(file.path(directory, "fit.rds"))
  stopifnot(!length(intersect(fit$recipe$training_strains, evaluation$strain_id)))
  features <- y1011_features(evaluation, panel, fit$recipe)
  base <- drop(ep_predict(features, fit$beta))
  saved <- fread(file.path(directory, "predictions.csv"))
  error <- 0
  for (j in seq_along(fit$quantiles)) {
    value <- fit$quantiles[j]
    current <- saved[tau == value]
    stopifnot(identical(current$sample_id, evaluation$sample_id), identical(current$strain_id, evaluation$strain_id))
    prediction <- base + fit$calibration[[j]][evaluation$condition]
    error <- max(error, abs(current$y - evaluation$y), abs(current$prediction - prediction),
                 abs(current$loss - ep_loss(evaluation$y, prediction, value)))
  }
  stopifnot(is.finite(error), error < 1e-9)
  checks[[length(checks) + 1L]] <- data.table(method = request$method, run = basename(directory), maximum_error = error)
  losses <- saved[, .(loss = mean(loss)), by = .(strain_id, tau)]
  losses[, method := request$method]
  records[[length(records) + 1L]] <- losses
  condition_loss <- saved[, .(strains = uniqueN(strain_id), loss = mean(loss),
                              coverage = mean(y <= prediction)), by = .(tau, condition)]
  condition_loss[, method := request$method]
  condition_records[[length(condition_records) + 1L]] <- condition_loss
  if (request$method == "JOINT") {
    active <- which(abs(fit$beta) > 1e-10)
    selected <- data.table(variable = names(fit$beta)[active], coefficient_scaled = fit$beta[active],
                           coefficient_response_scale = fit$beta[active] * fit$recipe$y_scale,
                           tau = request$tau, role = ifelse(features$penalty_factor[active] == 0, "adjustment", "selected"))
    setting <- as.data.frame(fit$grid)[fit$choice, , drop = FALSE]
    selected[, `:=`(Lambda = setting$Lambda, size = setting$size)]
    selections[[length(selections) + 1L]] <- selected
  }
  rm(features, fit)
  gc()
}
losses <- rbindlist(records)
stopifnot(!anyDuplicated(losses[, .(strain_id, tau, method)]))
summary <- losses[, .(strains = .N, loss = mean(loss), se = sd(loss) / sqrt(.N)), by = .(method, tau)]
comparisons <- list()
set.seed(20261104L)
for (value in sort(unique(losses$tau))) {
  wide <- dcast(losses[tau == value], strain_id ~ method, value.var = "loss")
  stopifnot("JOINT" %in% names(wide), !anyNA(wide))
  n <- nrow(wide)
  weights <- rmultinom(10000L, n, rep(1 / n, n)) / n
  for (method_name in setdiff(names(wide), c("strain_id", "JOINT"))) {
    reference <- wide[[method_name]]
    difference <- reference - wide$JOINT
    bootstrap <- drop(crossprod(weights, difference)) / drop(crossprod(weights, reference))
    comparisons[[length(comparisons) + 1L]] <- data.table(tau = value, comparator = method_name,
      strains = n, joint_loss = mean(wide$JOINT), comparator_loss = mean(reference),
      relative_gain = mean(difference) / mean(reference), lower = quantile(bootstrap, 0.025),
      upper = quantile(bootstrap, 0.975), phase = "discovery", manuscript_adoption = FALSE)
  }
}
output <- "results/preflight/empirical_selection_20261101/yeast_1011_prediction"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
comparison <- rbindlist(comparisons)
fwrite(summary, file.path(output, "summary.csv"))
fwrite(rbindlist(checks), file.path(output, "artifact_checks.csv"))
fwrite(rbindlist(condition_records), file.path(output, "condition_summary.csv"))
fwrite(rbindlist(selections), file.path(output, "selected_coefficients.csv"))
fwrite(comparison, file.path(output, "comparisons.csv"))
write_json(list(stage = "discovery", completed_runs = nrow(queue), decision_quantile = config$decision_quantile,
                formal_validation_used = FALSE, manuscript_adoption = FALSE),
           file.path(output, "status.json"), auto_unbox = TRUE, pretty = TRUE)
print(comparison)
