source("research/submission_completion_20261027/load.R")
registry <- submission_registry("main")
args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) args[1L] else "results/scratch/submission_completion_20261027/variant_checks"
if (dir.exists(output)) stop("The variant verification directory already exists.")
dir.create(output, recursive = TRUE)
records <- list()
for (id in c("A01", "C03", "C14", "B16", "F01", "E07", "E08", "G05")) {
  config <- config_row_to_list_v2(registry[registry$experiment_id == id, ])
  config$methods <- config$methods[startsWith(config$methods, "JOINT-PROFILE")]
  result <- run_one_replication_v2(".", config, 1L, base_seed = 202610271L,
    final = TRUE, allow_missing_benchmarks = FALSE, run_id = "submission_variant_checks")
  saveRDS(result, file.path(output, paste0(id, ".rds")))
  validate_paired_failed_coordinate_accounting_v2(result$metrics, result$coordinates)
  if (any(!result$metrics$converged) || any(!is.finite(result$coordinates$estimated_se))) {
    print(result$metrics[c("method_id", "status", "failure_message")])
    stop("A joint procedure did not pass its numerical verification.")
  }
  if (nrow(result$coordinates) != length(config$methods) * config$target_coordinate_count) stop("Incomplete joint inference coordinates.")
  if (id == "F01") {
    if (result$screening$overlap_count != 0L) stop("Independent screening subjects overlap.")
    method <- "JOINT-PROFILE-SPLIT"
    beta <- map_beta_to_full_v2(result$method_diagnostics[[method]]$beta_hat,
      names(make_beta_v2(config$p, config$s, config$signal)))
    expected <- sqrt(sum((beta - make_beta_v2(config$p, config$s, config$signal))^2))
    observed <- result$metrics$l2_error[result$metrics$method_id == method]
    if (abs(expected - observed) > 1e-10) stop("Screened estimation error omits a registered coordinate.")
  }
  records[[length(records) + 1L]] <- data.frame(experiment_id = id, methods = nrow(result$metrics),
    coordinates = nrow(result$coordinates), maximum_kkt = max(result$metrics$kkt_residual, na.rm = TRUE),
    maximum_runtime = max(result$metrics$runtime_sec))
  cat("Verified", id, "\n")
}
write.csv(do.call(rbind, records), file.path(output, "checks.csv"), row.names = FALSE)
print(do.call(rbind, records))
