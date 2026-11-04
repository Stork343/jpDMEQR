library(data.table)
library(jsonlite)
library(yaml)
library(BEDMatrix)
source("research/arabidopsis_application_20261105/features.R")
source("research/empirical_selection_20261101/prediction_core.R")

config <- read_yaml("config/discovery_arabidopsis_20261105.yml")
queue <- fread("results/raw/application/arabidopsis_20261105/screening_20261105/progress.csv")
stopifnot(nrow(queue) == 5L, all(queue$status == "completed"))
panel <- arab_panel(config$input_directory)
evaluation <- fread(file.path(config$input_directory, "discovery.csv"), colClasses = c(accession_id = "character"))[exploration_role == "evaluation"]
stopifnot(all(evaluation$role == "discovery"))
records <- checks <- selected_effects <- settings <- list()
for (j in seq_len(nrow(queue))) {
  directory <- queue$directory[j]
  request <- read_json(file.path(directory, "request.json"), simplifyVector = TRUE)
  status <- read_json(file.path(directory, "status.json"), simplifyVector = TRUE)
  stopifnot(status$status == "completed", !request$formal_validation_used)
  fit <- readRDS(file.path(directory, "fit.rds"))
  stopifnot(!length(intersect(fit$recipe$training_accessions, evaluation$accession_id)))
  features <- arab_features(evaluation, panel, fit$recipe)
  base <- drop(ep_predict(features, fit$beta))
  saved <- fread(file.path(directory, "predictions.csv"), colClasses = c(accession_id = "character"))
  error <- 0
  for (k in seq_along(fit$quantiles)) {
    value <- fit$quantiles[k]
    current <- saved[tau == value]
    stopifnot(identical(current$sample_id, evaluation$sample_id), identical(current$accession_id, evaluation$accession_id))
    prediction <- base + fit$calibration[[k]][evaluation$condition]
    error <- max(error, abs(current$y - evaluation$y), abs(current$prediction - prediction),
                 abs(current$loss - ep_loss(evaluation$y, prediction, value)))
  }
  stopifnot(is.finite(error), error < 1e-9)
  checks[[j]] <- data.table(method = request$method, run = basename(directory), maximum_error = error)
  losses <- saved[, .(loss = mean(loss)), by = .(accession_id, tau)]
  losses[, method := request$method]
  records[[j]] <- losses
  setting <- as.data.frame(fit$grid)[fit$choice, , drop = FALSE]
  settings[[j]] <- cbind(data.table(method = request$method, tau = request$tau,
    runtime_sec = status$runtime_sec), as.data.table(setting))
  if (request$method == "JOINT") {
    d <- length(fit$recipe$low_columns)
    p <- length(fit$recipe$markers)
    main <- fit$beta[d + seq_len(p)]
    interaction <- fit$beta[d + p + seq_len(p)]
    active <- which(abs(main) > 1e-10 | abs(interaction) > 1e-10)
    selected_effects[[length(selected_effects) + 1L]] <- data.table(
      marker = fit$recipe$markers[active], tau = request$tau,
      FT10_days = 2 * (main[active] - interaction[active]) * fit$recipe$y_scale / fit$recipe$scale[active],
      FT16_days = 2 * (main[active] + interaction[active]) * fit$recipe$y_scale / fit$recipe$scale[active],
      temperature_difference_days = 4 * interaction[active] * fit$recipe$y_scale / fit$recipe$scale[active],
      contrast = "source_binary_state_1_minus_0", phase = "discovery")
  }
  rm(features, fit)
  gc()
}
losses <- rbindlist(records)
stopifnot(!anyDuplicated(losses[, .(accession_id, tau, method)]))
summary <- losses[, .(accessions = .N, loss = mean(loss), se = sd(loss) / sqrt(.N)), by = .(method, tau)]
set.seed(20261105L)
comparisons <- list()
for (value in sort(unique(losses$tau))) {
  wide <- dcast(losses[tau == value], accession_id ~ method, value.var = "loss")
  stopifnot("JOINT" %in% names(wide), !anyNA(wide), nrow(wide) == 170L)
  weights <- rmultinom(10000L, nrow(wide), rep(1 / nrow(wide), nrow(wide))) / nrow(wide)
  for (method_name in setdiff(names(wide), c("accession_id", "JOINT"))) {
    reference <- wide[[method_name]]
    difference <- reference - wide$JOINT
    bootstrap <- drop(crossprod(weights, difference)) / drop(crossprod(weights, reference))
    comparisons[[length(comparisons) + 1L]] <- data.table(tau = value, comparator = method_name,
      accessions = nrow(wide), joint_loss = mean(wide$JOINT), comparator_loss = mean(reference),
      relative_gain = mean(difference) / mean(reference), lower = quantile(bootstrap, 0.025),
      upper = quantile(bootstrap, 0.975), phase = "discovery", manuscript_adoption = FALSE)
  }
}
output <- "results/preflight/empirical_selection_20261101/arabidopsis_prediction"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
comparison <- rbindlist(comparisons)
fwrite(summary, file.path(output, "summary.csv"))
fwrite(rbindlist(checks), file.path(output, "artifact_checks.csv"))
fwrite(rbindlist(settings, fill = TRUE), file.path(output, "settings.csv"))
fwrite(rbindlist(selected_effects), file.path(output, "selected_effects.csv"))
fwrite(comparison, file.path(output, "comparisons.csv"))
write_json(list(stage = "discovery", completed_runs = nrow(queue), decision_quantile = config$decision_quantile,
                formal_validation_used = FALSE, manuscript_adoption = FALSE),
           file.path(output, "status.json"), auto_unbox = TRUE, pretty = TRUE)
print(comparison)
