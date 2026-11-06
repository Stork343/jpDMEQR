library(yaml)
library(jsonlite)
source("research/clinical_quantile_20261105/engine.R")
source("research/clinical_quantile_20261105/orthogonal.R")
args <- commandArgs(trailingOnly = TRUE)
config <- read_yaml(if (length(args)) args[1L] else "config/clinical_quantile_smooth_20261106.yml")
checks <- list()
for (path in unlist(config$source_studies)) {
  study <- read_yaml(path)
  root <- file.path(config$output_root, paste0("orthogonal_", study$study_id))
  if (!dir.exists(root)) next
  status <- read_json(file.path(root, "status.json"), simplifyVector = TRUE)
  stopifnot(status$status == "completed")
  panel <- cp_panel(study$input_directory)
  reconstructed <- list()
  for (fold in seq_len(study$outer_folds)) {
    rows <- fread(file.path(study$output_root, paste0("outer_", fold), "evaluation.csv"))
    stopifnot(all(rows$role == "discovery_evaluation"))
    central <- readRDS(file.path(config$output_root, study$study_id, paste0("outer_", fold),
      sprintf("tau%03d", as.integer(100 * status$tau)), "fit.rds"))
    stopifnot(!length(intersect(central$recipe$training_strains, rows$strain_id)))
    features <- cq_features(rows, panel, central$recipe, central$columns)
    nuisance <- readRDS(file.path(root, paste0("outer_", fold), "nuisance.rds"))
    density <- cq_density(features, nuisance$density_fits[[1L]]$beta,
      nuisance$density_fits[[2L]]$beta, nuisance$delta, config)
    for (method in c("independence", "exchangeable")) {
      saved <- readRDS(file.path(root, paste0("outer_", fold), paste0(method, ".rds")))
      result <- cq_influence(features, central$beta, density$value, saved$projection)
      stopifnot(identical(result$patient, saved$result$patient))
      error <- max(abs(result$influence - saved$result$influence), abs(result$initial - saved$result$initial))
      stopifnot(error < 1e-10)
      checks[[length(checks) + 1L]] <- data.table(study_id = study$study_id, fold, method,
        patients = length(result$patient), reconstruction_error = error)
      reconstructed[[length(reconstructed) + 1L]] <- data.table(fold, method,
        strain_id = result$patient, initial = result$initial, influence = result$influence)
    }
  }
  all <- rbindlist(reconstructed)
  all[, centered := influence - mean(influence), by = .(fold, method)]
  summary <- all[, .(patients = .N, estimate = mean(initial + influence),
    standard_error = sqrt(sum(centered^2)) / .N), by = method]
  stored <- fread(file.path(config$report_directory, paste0(study$study_id, "_orthogonal_summary.csv")))
  comparison <- merge(summary, stored, by = "method", suffixes = c("_reconstructed", "_stored"))
  stopifnot(nrow(comparison) == 2L,
    max(abs(comparison$estimate_reconstructed - comparison$estimate_stored)) < 1e-10,
    max(abs(comparison$standard_error_reconstructed - comparison$standard_error_stored)) < 1e-10)
}
expected <- if (length(args)) 2L * sum(vapply(unlist(config$source_studies),
  function(path) read_yaml(path)$outer_folds, integer(1))) else 20L
stopifnot(length(checks) == expected)
fwrite(rbindlist(checks), file.path(config$report_directory, "influence_reconstruction_checks.csv"))
write_json(list(status = "passed", fitted_directions = length(checks),
  maximum_influence_error = max(rbindlist(checks)$reconstruction_error),
  studies = unique(rbindlist(checks)$study_id)), file.path(config$report_directory, "influence_reconstruction.json"),
  auto_unbox = TRUE, pretty = TRUE)
print(rbindlist(checks)[, .(fitted_directions = .N, maximum_error = max(reconstruction_error)), by = study_id])
