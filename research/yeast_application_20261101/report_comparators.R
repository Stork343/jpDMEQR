library(data.table)
library(jsonlite)
library(yaml)
source("research/yeast_application_20261101/prediction.R")
source("research/yeast_application_20261101/comparators.R")

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
config <- read_yaml(args[1L])
trait_name <- args[2L]
panel <- readRDS(config$input_panel)
data <- readRDS(config$input_growth)
stopifnot(all(data$role == "discovery"))
evaluation <- data[trait == trait_name & exploration_role == "evaluation"]
stopifnot(nrow(evaluation) > 0L, uniqueN(evaluation$trait) == 1L)
directories <- list.dirs(file.path(config$output_root, trait_name), recursive = FALSE)
directories <- directories[file.exists(file.path(directories, "status.json"))]
all_losses <- list()
checks <- list()
for (directory in directories) {
  request <- read_json(file.path(directory, "request.json"), simplifyVector = TRUE)
  if (request$method == "ENET" && basename(directory) != paste0("ENET_", config$mean_run_id)) next
  status <- read_json(file.path(directory, "status.json"), simplifyVector = TRUE)
  stopifnot(status$status == "completed", !request$formal_validation_used)
  fit <- readRDS(file.path(directory, "fit.rds"))
  saved <- fread(file.path(directory, "predictions.csv"))
  stopifnot(!length(intersect(fit$recipe$training_strains, evaluation$strain_id)))
  features <- yeast_features(evaluation, panel, fit$recipe)
  quantiles <- sort(unique(saved$tau))
  if (request$method == "FOREST") {
    library(grf)
    base <- matrix(NA_real_, nrow(evaluation), length(quantiles))
    for (g in unique(fit$selected)) {
      current <- yeast_forest_predict(fit$models[[as.character(g)]], features)
      for (k in which(fit$selected == g)) base[, k] <- current[, k]
    }
  } else base <- matrix(yeast_raw_prediction(features, fit$beta), nrow(evaluation), length(quantiles))
  error <- 0
  for (k in seq_along(quantiles)) {
    tau_value <- quantiles[k]
    calibration <- if (request$method %in% c("JOINT", "POOLED")) fit$calibration else fit$calibration[[k]]
    predicted <- base[, k] + calibration[as.character(evaluation$replicate)]
    observed <- saved[tau == tau_value]
    stopifnot(identical(observed$strain_id, evaluation$strain_id),
              identical(as.integer(observed$replicate), as.integer(evaluation$replicate)))
    error <- max(error, abs(observed$y - evaluation$y), abs(predicted - observed$prediction),
                 abs(yeast_loss(evaluation$y, predicted, tau_value) - observed$loss))
  }
  stopifnot(is.finite(error), error < 1e-9)
  losses <- saved[, .(loss = mean(loss)), by = .(strain_id, tau)]
  losses[, method := request$method]
  all_losses[[length(all_losses) + 1L]] <- losses
  checks[[length(checks) + 1L]] <- data.table(trait = trait_name, method = request$method,
                                            run = basename(directory), maximum_error = error)
  rm(fit, features)
  gc()
}
losses <- rbindlist(all_losses)
stopifnot(!anyDuplicated(losses[, .(strain_id, tau, method)]))
summary <- losses[, .(strains = .N, loss = mean(loss), se = sd(loss) / sqrt(.N)), by = .(method, tau)]
comparisons <- list()
set.seed(20261202L)
for (tau_value in sort(unique(losses$tau))) {
  current <- dcast(losses[tau == tau_value], strain_id ~ method, value.var = "loss")
  if (!("JOINT" %in% names(current))) next
  stopifnot(!anyNA(current))
  n <- nrow(current)
  weights <- rmultinom(10000L, n, rep(1 / n, n)) / n
  for (method_name in setdiff(names(current), c("strain_id", "JOINT"))) {
    reference <- current[[method_name]]
    difference <- reference - current$JOINT
    bootstrap_gain <- drop(crossprod(weights, difference)) / drop(crossprod(weights, reference))
    comparisons[[length(comparisons) + 1L]] <- data.table(trait = trait_name, tau = tau_value,
      comparator = method_name, strains = n, joint_loss = mean(current$JOINT), comparator_loss = mean(reference),
      relative_gain = mean(difference) / mean(reference), lower = quantile(bootstrap_gain, 0.025),
      upper = quantile(bootstrap_gain, 0.975), phase = "discovery", manuscript_adoption = FALSE)
  }
}
output <- file.path("results/preflight/empirical_selection_20261101/yeast_comparators", basename(dirname(config$input_panel)), trait_name)
dir.create(output, recursive = TRUE, showWarnings = FALSE)
fwrite(summary, file.path(output, "summary.csv"))
fwrite(rbindlist(checks), file.path(output, "artifact_checks.csv"))
if (length(comparisons)) {
  comparison <- rbindlist(comparisons)
  fwrite(comparison, file.path(output, "comparisons.csv"))
  print(comparison)
}
write_json(list(stage = "discovery", completed_runs = length(all_losses), mean_run_id = config$mean_run_id,
                formal_validation_used = FALSE, manuscript_adoption = FALSE),
           file.path(output, "status.json"), auto_unbox = TRUE, pretty = TRUE)
