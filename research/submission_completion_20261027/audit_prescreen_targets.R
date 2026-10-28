library(data.table)
source("scripts/00_source_v2.R")
source_v2_module(".", "profile_v2", .GlobalEnv)
source_v2_module(".", "simulation_v2", .GlobalEnv)
root <- "results/raw/simulation/submission_20261027_main_joint"
if (!file.exists(file.path(root, "completed.rds"))) stop("The complete proposed-method study is required.")
request <- readRDS(file.path(root, "request.rds"))
jobs <- request$jobs[request$jobs$experiment_id %in% c("F05", "F06"), ]
records <- list()
for (i in seq_len(nrow(jobs))) {
  job <- jobs[i, ]
  path <- file.path(root, "replications", job$experiment_id, sprintf("replicate_%05d/result.rds", job$replicate))
  result <- readRDS(path)
  configuration <- result$config
  metric <- result$metrics[result$metrics$method_id == "JOINT-PROFILE", ]
  if (!isTRUE(metric$converged) || configuration$design_type != "ar1" ||
      configuration$x_b_corr != 0 || configuration$informative_size || configuration$heteroskedastic)
    stop("The Gaussian projection audit requires its declared generating distribution.")
  beta <- make_beta_v2(configuration$p, configuration$s, configuration$signal)
  diagnostic <- result$method_diagnostics[["JOINT-PROFILE"]]
  working_names <- names(diagnostic$beta_hat)
  if (!identical(working_names, result$screening$selected_names)) stop("The saved prescreening coordinate order differs.")
  fit <- diagnostic$submission$fit
  for (j in seq_len(nrow(fit$table))) for (fold in 1:2) {
    training <- fit$training[[fold]]
    support <- match(working_names[training$graphs[[j]]$support], names(beta))
    target <- match(fit$table$coordinate[j], names(beta))
    if (anyNA(support) || support[1L] != target) stop("A selected target coordinate cannot be reconstructed.")
    covariance <- configuration$rho_x^abs(outer(support, support, "-"))
    active <- which(beta != 0)
    projection <- drop(solve(covariance, configuration$rho_x^abs(outer(support, active, "-")) %*% beta[active]))
    records[[length(records) + 1L]] <- data.table(experiment_id = job$experiment_id, replicate = job$replicate,
      coordinate = fit$table$coordinate[j], coordinate_type = if (target %in% active) "active" else "null", fold = fold,
      estimation_subjects = length(training$validation_subjects), available_predictors = length(working_names),
      selected_predictors = length(support), full_target = beta[target], selected_target = projection[1L],
      displacement = projection[1L] - beta[target])
  }
}
records <- rbindlist(records)
paired <- records[, .(displacement = weighted.mean(displacement, estimation_subjects),
  maximum_fold_displacement = max(abs(displacement))), by = .(experiment_id, replicate, coordinate, coordinate_type)]
summary <- paired[, {
  shifted <- abs(displacement) > 1e-8
  interval <- binom.test(sum(shifted), .N)$conf.int
  list(replications = .N, shifted = sum(shifted), fraction_shifted = mean(shifted),
    fraction_mcse = sd(as.numeric(shifted)) / sqrt(.N), fraction_lower = interval[1L], fraction_upper = interval[2L],
    mean_displacement = mean(displacement), displacement_mcse = sd(displacement) / sqrt(.N),
    maximum_absolute_displacement = max(abs(displacement)))
}, by = .(experiment_id, coordinate, coordinate_type)]
output <- "results/tables/submission_20261028/prescreen_targets"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
fwrite(records, file.path(output, "fold_target_displacements.csv"))
fwrite(paired, file.path(output, "combined_target_displacements.csv"))
fwrite(summary, file.path(output, "coordinate_summary.csv"))
jsonlite::write_json(list(source = root, settings = c("F05", "F06"), datasets = nrow(jobs),
  estimator_fits_recomputed = FALSE, target_formula = "exact_Gaussian_conditional_projection",
  full_target_coverage_unchanged = TRUE, generated = as.character(Sys.time())),
  file.path(output, "manifest.json"), auto_unbox = TRUE, pretty = TRUE)
print(summary)
