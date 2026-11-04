library(data.table)
library(jsonlite)
library(yaml)
library(BEDMatrix)
source("research/rice_baap_application_20261104/features.R")
source("research/empirical_selection_20261101/prediction_core.R")

rice_saved_prediction <- function(rows, panel, fit) {
  recipe <- fit$recipe
  p <- length(recipe$markers)
  d <- length(recipe$low_columns)
  k <- ncol(recipe$loadings)
  stopifnot(length(fit$beta) == d + 2L * p)
  ids <- sort(unique(rows$accession_id))
  high <- fit$beta[d + seq_len(p)]
  interaction <- fit$beta[d + p + seq_len(p)] / recipe$interaction_scale
  projection <- cbind(recipe$loadings, marker_effect = high, water_effect = interaction)
  values <- matrix(0, length(ids), k + 2L, dimnames = list(ids, colnames(projection)))
  bed <- panel$beds[[recipe$panel_key]]
  for (start in seq.int(1L, p, by = 4096L)) {
    j <- start:min(start + 4095L, p)
    G <- bed[ids, recipe$markers[j], drop = FALSE]
    G <- sweep(sweep(G, 2L, recipe$center[j], "-"), 2L, recipe$scale[j], "/")
    values <- values + G %*% projection[j, , drop = FALSE]
  }
  low <- rice_low(rows, values[, seq_len(k), drop = FALSE])[, recipe$low_columns, drop = FALSE]
  normalized <- drop(low %*% fit$beta[seq_len(d)]) +
    values[rows$accession_id, "marker_effect"] +
    (as.numeric(rows$irrigation == "AWD") - 0.5) * values[rows$accession_id, "water_effect"] -
    sum(recipe$interaction_center * interaction)
  normalized * recipe$y_scale + recipe$y_center
}

config <- read_yaml("config/discovery_rice_baap_20261104.yml")
panel <- rice_panel(config$input_directory)
evaluation <- fread(file.path(config$input_directory, "discovery.csv"))[exploration_role == "evaluation"]
stopifnot(all(evaluation$role == "discovery"))
directories <- list.dirs(file.path(config$output_root, "Cadmium"), recursive = FALSE)
directories <- directories[file.exists(file.path(directories, "status.json"))]
records <- checks <- list()
for (directory in directories) {
  request <- read_json(file.path(directory, "request.json"), simplifyVector = TRUE)
  if (request$method == "ENET" && basename(directory) != paste0("ENET_", config$mean_run_id)) next
  if (request$method == "JOINT" && !startsWith(basename(directory), paste0("JOINT_", config$joint_run_prefix, "_"))) next
  status <- read_json(file.path(directory, "status.json"), simplifyVector = TRUE)
  stopifnot(status$status == "completed", !request$formal_validation_used)
  fit <- readRDS(file.path(directory, "fit.rds"))
  stopifnot(!length(intersect(fit$recipe$training_accessions, evaluation$accession_id)))
  base <- rice_saved_prediction(evaluation, panel, fit)
  saved <- fread(file.path(directory, "predictions.csv"))
  error <- 0
  for (j in seq_along(fit$quantiles)) {
    value <- fit$quantiles[j]
    current <- saved[tau == value]
    stopifnot(identical(current$sample_id, evaluation$sample_id), identical(current$accession_id, evaluation$accession_id))
    prediction <- base + fit$calibration[[j]][evaluation$environment]
    error <- max(error, abs(current$y - evaluation$y), abs(current$prediction - prediction),
                 abs(current$loss - ep_loss(evaluation$y, prediction, value)))
  }
  stopifnot(is.finite(error), error < 1e-9)
  checks[[length(checks) + 1L]] <- data.table(method = request$method, run = basename(directory), maximum_error = error)
  losses <- saved[, .(loss = mean(loss)), by = .(accession_id, tau)]
  losses[, method := request$method]
  records[[length(records) + 1L]] <- losses
  rm(fit)
  gc()
}
stopifnot(length(records) > 0L)
losses <- rbindlist(records)
stopifnot(!anyDuplicated(losses[, .(accession_id, tau, method)]))
summary <- losses[, .(accessions = .N, loss = mean(loss), se = sd(loss) / sqrt(.N)), by = .(method, tau)]
comparisons <- list()
set.seed(20261104L)
for (value in sort(unique(losses$tau))) {
  wide <- dcast(losses[tau == value], accession_id ~ method, value.var = "loss")
  if (!("JOINT" %in% names(wide))) next
  stopifnot(!anyNA(wide))
  n <- nrow(wide)
  weights <- rmultinom(10000L, n, rep(1 / n, n)) / n
  for (method_name in setdiff(names(wide), c("accession_id", "JOINT"))) {
    reference <- wide[[method_name]]
    difference <- reference - wide$JOINT
    bootstrap <- drop(crossprod(weights, difference)) / drop(crossprod(weights, reference))
    comparisons[[length(comparisons) + 1L]] <- data.table(tau = value, comparator = method_name,
      accessions = n, joint_loss = mean(wide$JOINT), comparator_loss = mean(reference),
      relative_gain = mean(difference) / mean(reference), lower = quantile(bootstrap, 0.025),
      upper = quantile(bootstrap, 0.975), phase = "discovery", manuscript_adoption = FALSE)
  }
}
output <- "results/preflight/empirical_selection_20261101/rice_baap_prediction"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
fwrite(summary, file.path(output, "summary.csv"))
fwrite(rbindlist(checks), file.path(output, "artifact_checks.csv"))
if (length(comparisons)) {
  comparison <- rbindlist(comparisons)
  fwrite(comparison, file.path(output, "comparisons.csv"))
  print(comparison)
}
write_json(list(stage = "discovery", completed_runs = length(records), decision_quantile = config$decision_quantile,
                formal_validation_used = FALSE, manuscript_adoption = FALSE),
           file.path(output, "status.json"), auto_unbox = TRUE, pretty = TRUE)
