library(data.table)
library(jsonlite)
library(yaml)
source("research/yeast_application_20261101/prediction.R")

root <- "results/raw/application/yeast_20261101/screening_20261101"
jobs <- fread(file.path(root, "progress.csv"))
config <- read_yaml("config/discovery_yeast_20261101.yml")
jobs[method == "ENET", `:=`(run_id = config$mean_run_id,
  directory = file.path(config$output_root, trait, paste0("ENET_", config$mean_run_id)))]
stopifnot(all(file.exists(file.path(jobs$directory, "status.json"))))
complete <- jobs[, .(finished = all(status == "completed")), by = trait][finished == TRUE, trait]
requested <- commandArgs(trailingOnly = TRUE)
if (length(requested)) complete <- intersect(complete, requested)
panel <- readRDS("data/derived/yeast_20261101/panel.rds")
data <- readRDS("data/derived/yeast_20261101/growth_discovery.rds")
output_dir <- "results/preflight/empirical_selection_20261101/yeast_screening"
if (length(requested)) output_dir <- file.path(output_dir, "selected_cases", paste(requested, collapse = "_"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
summaries <- list()
checks <- list()
for (trait_name in complete) {
  current_jobs <- jobs[trait == trait_name]
  evaluation <- data[trait == trait_name & exploration_role == "evaluation"]
  all_losses <- list()
  for (j in seq_len(nrow(current_jobs))) {
    directory <- current_jobs$directory[j]
    fit <- readRDS(file.path(directory, "fit.rds"))
    stopifnot(!length(intersect(fit$recipe$training_strains, evaluation$strain_id)))
    features <- yeast_features(evaluation, panel, fit$recipe)
    base_prediction <- drop(yeast_raw_prediction(features, fit$beta))
    saved <- fread(file.path(directory, "predictions.csv"))
    quantiles <- sort(unique(saved$tau))
    maximum_error <- 0
    for (k in seq_along(quantiles)) {
      tau_value <- quantiles[k]
      calibration <- if (current_jobs$method[j] == "JOINT") fit$calibration else fit$calibration[[k]]
      predicted <- base_prediction + calibration[as.character(evaluation$replicate)]
      observed <- saved[tau == tau_value]
      stopifnot(identical(observed$strain_id, evaluation$strain_id),
                identical(as.integer(observed$replicate), as.integer(evaluation$replicate)))
      error <- max(abs(observed$y - evaluation$y), abs(predicted - observed$prediction),
                   abs(yeast_loss(evaluation$y, predicted, tau_value) - observed$loss))
      maximum_error <- max(maximum_error, error)
    }
    stopifnot(maximum_error < 1e-9)
    report <- list(status = "prediction_artifacts_verified", maximum_error = maximum_error,
                   evaluation_strains = uniqueN(evaluation$strain_id), training_overlap = 0L,
                   formal_validation_used = FALSE)
    write_json(report, file.path(directory, "artifact_checks.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
    checks[[length(checks) + 1L]] <- data.table(trait = trait_name, method = current_jobs$method[j],
                                               run_id = current_jobs$run_id[j], maximum_error = maximum_error)
    losses <- saved[, .(loss = mean(loss)), by = .(strain_id, tau)]
    losses[, method := current_jobs$method[j]]
    all_losses[[j]] <- losses
    rm(features, fit)
    gc()
  }
  wide <- dcast(rbindlist(all_losses), strain_id + tau ~ method, value.var = "loss")
  stopifnot(!anyNA(wide), !anyDuplicated(wide[, .(strain_id, tau)]))
  set.seed(20261130L + match(trait_name, sort(unique(jobs$trait))))
  for (tau_value in sort(unique(wide$tau))) {
    values <- wide[tau == tau_value]
    n <- nrow(values)
    weights <- rmultinom(10000L, n, rep(1 / n, n)) / n
    difference <- values$ENET - values$JOINT
    bootstrap_gain <- drop(crossprod(weights, difference)) / drop(crossprod(weights, values$ENET))
    summaries[[length(summaries) + 1L]] <- data.table(
      trait = trait_name, tau = tau_value, strains = n,
      joint_loss = mean(values$JOINT), enet_loss = mean(values$ENET),
      relative_gain = mean(difference) / mean(values$ENET),
      lower = quantile(bootstrap_gain, 0.025), upper = quantile(bootstrap_gain, 0.975),
      phase = "discovery", manuscript_adoption = FALSE)
  }
}
if (length(summaries)) {
  summary <- rbindlist(summaries)
  fwrite(summary, file.path(output_dir, "comparisons.csv"))
  fwrite(rbindlist(checks), file.path(output_dir, "artifact_checks.csv"))
  fwrite(summary[order(-relative_gain)], file.path(output_dir, "ranked_comparisons.csv"))
  print(summary[tau == 0.25][order(-relative_gain)][1:min(10L, .N)])
}
write_json(list(complete_conditions = length(complete), total_conditions = uniqueN(jobs$trait),
                formal_validation_used = FALSE, manuscript_adoption = FALSE,
                updated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)),
           file.path(output_dir, "status.json"), auto_unbox = TRUE, pretty = TRUE)
