library(data.table)
library(yaml)
library(jsonlite)
config <- read_yaml("config/application_soynam_20261031.yml")
input <- readRDS(file.path(config$results_directory, "input.rds"))
registry <- fread(file.path(config$results_directory, "registry.csv"))
inference_root <- file.path(config$results_directory, "inference", config$execution_version)
prediction_root <- file.path(config$results_directory, "prediction", config$execution_version)
records <- list()
for (case_id in registry$case_id) {
  directory <- file.path(inference_root, case_id)
  backup <- file.path(config$backup_directory, "inference", config$execution_version, case_id)
  required <- c("request.rds", "prepared.rds", "fit.rds", "inference.csv", "diagnostics.csv", "completion.json")
  present <- file.exists(file.path(directory, required)) & file.exists(file.path(backup, required))
  if (!all(present)) stop("An inference case or its backup is incomplete: ", case_id)
  if (any(file.info(file.path(directory, required))$size != file.info(file.path(backup, required))$size))
    stop("An inference backup has a different length: ", case_id)
  current <- fread(file.path(directory, "inference.csv"))
  diagnostics <- fread(file.path(directory, "diagnostics.csv"), integer64 = "double")
  if (nrow(current) != 2L || nrow(diagnostics) != 8L || any(diagnostics$kkt > 1e-8) ||
      any(diagnostics$nuisance_gradient > 1e-8) || any(diagnostics$basis_minimum <= 0))
    stop("An inference result fails its declared checks: ", case_id)
  records[[length(records) + 1L]] <- data.table(case_id, kind = "inference", status = "complete", backup = TRUE,
    n = unique(current$n), systems = nrow(diagnostics), maximum_kkt = max(diagnostics$kkt))
}
for (method in unlist(config$prediction$methods)) for (tau in unlist(config$quantiles)) {
  case_id <- paste0(method, "_q", sprintf("%02d", round(100 * tau)))
  directory <- file.path(prediction_root, case_id)
  backup <- file.path(config$backup_directory, "prediction", config$execution_version, case_id)
  required <- c("request.rds", "model.rds", "predictions.csv", "metrics.csv", "completion.json")
  if (!all(file.exists(file.path(directory, required)) & file.exists(file.path(backup, required))))
    stop("A prediction case or its backup is incomplete: ", case_id)
  if (any(file.info(file.path(directory, required))$size != file.info(file.path(backup, required))$size))
    stop("A prediction backup has a different length: ", case_id)
  request <- readRDS(file.path(directory, "request.rds"))
  model <- readRDS(file.path(directory, "model.rds"))
  data <- fread(file.path(directory, "predictions.csv"))
  expected <- input$line_folds[prediction_fold == config$prediction$heldout_fold, cluster_id]
  if (!setequal(data$cluster_id, expected) || length(intersect(request$training_lines, request$test_lines)))
    stop("The prediction sample allocation differs: ", case_id)
  if (any(!is.finite(data$prediction)) || length(model$oof) != 2L * length(request$training_lines) ||
      any(!is.finite(model$oof))) stop("The prediction output is incomplete: ", case_id)
  recalculated <- (data$observed - data$prediction) * (tau - (data$observed < data$prediction))
  if (max(abs(recalculated - data$loss)) > 1e-10) stop("The prediction loss differs from its recorded values.")
  if (method == "pooled_quantile_lasso") {
    paths <- c(paste0("training_path_", 1:3, ".rds"), "training_path_final.rds")
    if (!all(file.exists(file.path(directory, paths)) & file.exists(file.path(backup, paths))))
      stop("A pooled training path or its backup is missing.")
    for (path in paths) {
      fitted <- readRDS(file.path(directory, path))
      if (length(fitted$lambda) != config$prediction$pooled_nlambda ||
          max(abs(fitted$diagnostics$primal_dual_gap)) > config$prediction$pooled_gap_tolerance ||
          max(fitted$diagnostics$dual_violation) > config$prediction$pooled_dual_tolerance ||
          max(fitted$diagnostics$mapping_error) > 1e-9)
        stop("A pooled training path fails its complete certificate.")
    }
  }
  records[[length(records) + 1L]] <- data.table(case_id, kind = "prediction", status = "complete", backup = TRUE,
    n = uniqueN(data$cluster_id), systems = NA_integer_, maximum_kkt = NA_real_)
}
pooled <- fread("results/preflight/soynam_20261031/pooled_lp_checks.csv")
pooled_path <- fread("results/preflight/soynam_20261031/pooled_lp_path_checks.csv")
paired <- fread("results/preflight/soynam_20261031/pair_qp_checks.csv")
if (nrow(pooled) != 9L || max(pooled$objective_difference) > 2e-6 ||
    max(pooled$mapping_error) > 1e-9) stop("The pooled implementation checks do not pass.")
fwrite(rbindlist(records), "results/preflight/soynam_20261031/delivery_checks.csv")
manifest <- list(status = "complete", inference_cases = nrow(registry), prediction_cases = 9L,
  original_coordinate_systems = sum(rbindlist(records)$systems, na.rm = TRUE),
  all_backups_present = TRUE, source_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  pooled_solver = "HiGHS 1.15.1", pooled_check_cases = nrow(pooled),
  pooled_maximum_objective_difference = max(pooled$objective_difference),
  pooled_full_path_checks = pooled_path$path_points,
  paired_solver_check_cases = nrow(paired), checked = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE))
write_json(manifest, "results/preflight/soynam_20261031/delivery_manifest.json", auto_unbox = TRUE, pretty = TRUE)
writeLines(trimws(capture.output(sessionInfo()), which = "right"), "results/preflight/soynam_20261031/session.txt")
print(manifest)
