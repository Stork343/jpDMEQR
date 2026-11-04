library(data.table)
library(jsonlite)
library(yaml)
source("research/yeast_application_20261101/prediction.R")
source("research/drugmatrix_application_20261101/prediction.R")

config <- read_yaml("config/discovery_drugmatrix_20261101.yml")
panel <- readRDS(config$input_panel)
rows <- readRDS(config$input_growth)
evaluation <- rows[exploration_role == "evaluation"]
stopifnot(all(evaluation$role == "discovery"))
directories <- list.dirs(file.path(config$output_root, "ALT"), recursive = FALSE)
directories <- directories[file.exists(file.path(directories, "status.json"))]
records <- list()
checks <- list()
for (directory in directories) {
  request <- read_json(file.path(directory, "request.json"), simplifyVector = TRUE)
  if (request$method == "ENET" && basename(directory) != paste0("ENET_", config$mean_run_id)) next
  if (request$method == "JOINT" && !startsWith(basename(directory), paste0("JOINT_", config$joint_run_prefix, "_"))) next
  status <- read_json(file.path(directory, "status.json"), simplifyVector = TRUE)
  stopifnot(status$status == "completed", !request$formal_validation_used)
  fit <- readRDS(file.path(directory, "fit.rds"))
  stopifnot(!length(intersect(fit$recipe$training_compounds, evaluation$compound_id)))
  features <- dm_features(evaluation, panel, fit$recipe)
  base <- drop(dm_predict(features, fit$beta))
  saved <- fread(file.path(directory, "predictions.csv"))
  maximum_error <- 0
  for (j in seq_along(fit$quantiles)) {
    value <- fit$quantiles[j]
    current <- saved[tau == value]
    stopifnot(identical(current$sample_id, evaluation$sample_id),
              identical(current$compound_id, evaluation$compound_id))
    predicted <- base + fit$calibration[j]
    maximum_error <- max(maximum_error, abs(current$y - evaluation$y),
      abs(current$prediction - predicted), abs(current$loss - yeast_loss(evaluation$y, predicted, value)))
  }
  stopifnot(is.finite(maximum_error), maximum_error < 1e-9)
  checks[[length(checks) + 1L]] <- data.table(method = request$method, run = basename(directory),
    maximum_error = maximum_error, training_overlap = 0L)
  loss <- saved[, .(loss = mean(loss)), by = .(compound_id, tau)]
  loss[, method := request$method]
  records[[length(records) + 1L]] <- loss
}
stopifnot(length(records) > 0L)
losses <- rbindlist(records)
stopifnot(!anyDuplicated(losses[, .(compound_id, tau, method)]))
summary <- losses[, .(compounds = .N, loss = mean(loss), se = sd(loss) / sqrt(.N)), by = .(method, tau)]
comparisons <- list()
set.seed(20261205L)
for (value in sort(unique(losses$tau))) {
  wide <- dcast(losses[tau == value], compound_id ~ method, value.var = "loss")
  if (!("JOINT" %in% names(wide))) next
  stopifnot(!anyNA(wide))
  n <- nrow(wide)
  weights <- rmultinom(10000L, n, rep(1 / n, n)) / n
  for (method_name in setdiff(names(wide), c("compound_id", "JOINT"))) {
    reference <- wide[[method_name]]
    difference <- reference - wide$JOINT
    boot <- drop(crossprod(weights, difference)) / drop(crossprod(weights, reference))
    comparisons[[length(comparisons) + 1L]] <- data.table(tau = value, comparator = method_name,
      compounds = n, joint_loss = mean(wide$JOINT), comparator_loss = mean(reference),
      relative_gain = mean(difference) / mean(reference), lower = quantile(boot, 0.025),
      upper = quantile(boot, 0.975), phase = "discovery", manuscript_adoption = FALSE)
  }
}
output <- "results/preflight/empirical_selection_20261101/drugmatrix_prediction"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
fwrite(summary, file.path(output, "summary.csv"))
fwrite(rbindlist(checks), file.path(output, "artifact_checks.csv"))
if (length(comparisons)) {
  table <- rbindlist(comparisons)
  fwrite(table, file.path(output, "comparisons.csv"))
  print(table)
}
write_json(list(stage = "discovery", completed_runs = length(records), mean_run_id = config$mean_run_id,
                decision_quantile = config$decision_quantile,
                formal_validation_used = FALSE, manuscript_adoption = FALSE),
           file.path(output, "status.json"), auto_unbox = TRUE, pretty = TRUE)
