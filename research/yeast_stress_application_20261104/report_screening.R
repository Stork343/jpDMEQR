library(data.table)
library(jsonlite)
library(yaml)
library(BEDMatrix)
source("research/yeast_1011_application_20261104/features.R")
source("research/yeast_stress_application_20261104/features.R")
source("research/empirical_selection_20261101/prediction_core.R")

config <- read_yaml("config/discovery_yeast_stress_20261104.yml")
queue <- fread("results/raw/application/yeast_stress_20261104/screening_20261104/progress.csv")
stopifnot(nrow(queue) == 12L, all(queue$status == "completed"))
panel <- y1011_panel(config$genomic_directory)
records <- checks <- selections <- list()
for (j in seq_len(nrow(queue))) {
  directory <- queue$directory[j]
  assay_name <- queue$trait[j]
  evaluation <- fread(file.path(config$input_directory, assay_name, "discovery.csv"))[exploration_role == "evaluation"]
  request <- read_json(file.path(directory, "request.json"), simplifyVector = TRUE)
  status <- read_json(file.path(directory, "status.json"), simplifyVector = TRUE)
  stopifnot(status$status == "completed", !request$formal_validation_used, all(evaluation$role == "discovery"))
  fit <- readRDS(file.path(directory, "fit.rds"))
  stopifnot(!length(intersect(fit$recipe$training_strains, evaluation$strain_id)))
  features <- stress_features(evaluation, panel, fit$recipe)
  base <- drop(ep_predict(features, fit$beta))
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
  checks[[j]] <- data.table(assay = assay_name, method = request$method, maximum_error = error)
  losses <- saved[tau == config$decision_quantile, .(loss = mean(loss)), by = strain_id]
  losses[, `:=`(assay = assay_name, method = request$method)]
  records[[j]] <- losses
  if (request$method == "JOINT") {
    active <- which(abs(fit$beta) > 1e-10)
    selections[[assay_name]] <- data.table(assay = assay_name, variable = names(fit$beta)[active],
      coefficient_scaled = fit$beta[active], coefficient_response_scale = fit$beta[active] * fit$recipe$y_scale,
      role = ifelse(features$penalty_factor[active] == 0, "adjustment", "selected"))
  }
  rm(features, fit)
  gc()
}
losses <- rbindlist(records)
set.seed(20261104L)
comparisons <- rbindlist(lapply(names(config$assays), function(assay_name) {
  wide <- dcast(losses[assay == assay_name], strain_id ~ method, value.var = "loss")
  stopifnot(!anyNA(wide), nrow(wide) == 148L)
  weights <- rmultinom(10000L, nrow(wide), rep(1 / nrow(wide), nrow(wide))) / nrow(wide)
  difference <- wide$ENET - wide$JOINT
  bootstrap <- drop(crossprod(weights, difference)) / drop(crossprod(weights, wide$ENET))
  data.table(assay = assay_name, tau = config$decision_quantile, strains = nrow(wide),
    joint_loss = mean(wide$JOINT), comparator_loss = mean(wide$ENET),
    relative_gain = mean(difference) / mean(wide$ENET), lower = quantile(bootstrap, 0.025),
    upper = quantile(bootstrap, 0.975), exact_comparison_required = mean(difference) > 0,
    manuscript_adoption = FALSE)
}))
output <- "results/preflight/empirical_selection_20261101/yeast_stress_prediction"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
fwrite(comparisons, file.path(output, "comparisons.csv"))
fwrite(rbindlist(checks), file.path(output, "artifact_checks.csv"))
fwrite(rbindlist(selections), file.path(output, "selected_coefficients.csv"))
write_json(list(stage = "discovery", completed_runs = nrow(queue),
                exact_comparisons_required = sum(comparisons$exact_comparison_required),
                formal_validation_used = FALSE, manuscript_adoption = FALSE),
           file.path(output, "status.json"), auto_unbox = TRUE, pretty = TRUE)
print(comparisons)
