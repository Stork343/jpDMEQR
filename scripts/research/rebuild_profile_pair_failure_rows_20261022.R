root <- normalizePath(".", mustWork = TRUE)
source(file.path(root, "scripts", "00_source_v2.R"))
for (module in c("profile_v2", "simulation_v2", "metrics_v2")) source_v2_module(root, module)
source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"))
run_id <- "integration_20261022_02"
input <- file.path(root, "results", "raw", "simulation", run_id)
out <- file.path(root, "results", "scratch", "profile_pair_failure_accounting_20261022")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
objects <- readRDS(file.path(input, "raw_objects.rds"))
keys <- if (identical(objects$format, "replication_file_index_v1")) objects$files$key else
  vapply(objects, function(object) sprintf("%s_%05d", object$config$experiment_id,
                                           object$config$replicate %||% object$replicate), character(1))
metrics <- utils::read.csv(file.path(input, "replication_metrics.csv"))
coordinates <- utils::read.csv(file.path(input, "coordinate_metrics.csv"))
added <- list()
for (key in keys) {
  object <- read_registry_replication_v2(objects, key, base_directory = input)
  failures <- object$metrics[object$metrics$status == "failed" &
    object$metrics$score_definition %in% "profile_bandwidth_pair_2_half_minus_full", , drop = FALSE]
  if (!nrow(failures)) next
  if (!identical(object$config$target_mode, "structural")) {
    stop("This reconstruction requires the saved structural target definition.")
  }
  beta <- make_beta_v2(object$config$p, object$config$s, object$config$signal)
  target <- list(beta = beta, target_type = "structural", target_mc_se = 0,
                 target_mc_se_by_coordinate = stats::setNames(rep(0, length(beta)), names(beta)))
  for (i in seq_len(nrow(failures))) {
    method <- failures$method_id[i]
    partial <- object$method_diagnostics[[method]]$partial_fit
    requested <- object$method_diagnostics[[method]]$requested_target_coordinates %||%
      partial$target_coordinates
    if (is.null(requested) || !length(requested)) stop("The saved failed fit has no requested coordinates.")
    added[[length(added) + 1L]] <- failed_paired_coordinate_rows_v2(
      run_id, object$config, failures$replicate[i], failures$seed[i], method,
      requested, target, names(beta), seq_len(object$config$s), failures$implementation_commit[i]
    )
  }
}
reconstructed <- do.call(rbind, added)
if (!nrow(reconstructed)) stop("No saved paired failures were reconstructed.")
if (!setequal(names(coordinates), names(reconstructed))) stop("Coordinate schema differs from saved raw outputs.")
complete <- rbind(coordinates, reconstructed[, names(coordinates)])
statistics <- lapply(seq_len(nrow(complete)), function(i) {
  coordinate_wald_statistics_v2(complete$beta_tilde[i], complete$estimated_se[i], complete$target_value[i])
})
complete$wald_z <- vapply(statistics, `[[`, numeric(1), "wald_z")
complete$p_value <- vapply(statistics, `[[`, numeric(1), "p_value")
validate_paired_failed_coordinate_accounting_v2(metrics, complete)
utils::write.csv(reconstructed, file.path(out, "reconstructed_failed_coordinates.csv"), row.names = FALSE)
utils::write.csv(complete, file.path(out, "coordinate_metrics_complete.csv"), row.names = FALSE)
utils::write.csv(data.frame(source_run = run_id, saved_coordinate_rows = nrow(coordinates),
                            failed_method_rows = sum(metrics$status == "failed" &
                              metrics$score_definition %in% "profile_bandwidth_pair_2_half_minus_full"),
                            added_failed_coordinates = nrow(reconstructed),
                            complete_coordinate_rows = nrow(complete)),
                 file.path(out, "accounting.csv"), row.names = FALSE)
print(utils::read.csv(file.path(out, "accounting.csv")))
