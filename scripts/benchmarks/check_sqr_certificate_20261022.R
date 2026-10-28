args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[1] else ".", mustWork = TRUE)
.libPaths(c(file.path(root, "results", "scratch", "profile_pair_package_20261022", "library"), .libPaths()))
library(jpDMEQR)
source(file.path(root, "scripts", "00_source_v2.R"))
for (module in c("profile_v2", "simulation_v2", "benchmark_adapters_v2")) source_v2_module(root, module, environment())
source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"))
directory <- file.path(root, "results", "scratch", "profile_pair_benchmark_20261022",
                        paste0("P01_certificate_", format(Sys.time(), "%Y%m%dT%H%M%OS6")))
dir.create(directory, recursive = TRUE)
reference_directory <- file.path(root, "results", "scratch", "profile_pair_benchmark_20261022",
                                  "registered_sqr_P01_20261022T225500.901488")
request <- readRDS(file.path(reference_directory, "request.rds"))
reference <- readRDS(file.path(reference_directory, "result.rds"))
started <- proc.time()[[3L]]
precision <- solve_clime_matrix_iid_v2(reference$fit_object$hessian,
                                       reference$inference_object$precision$gamma_grid_requested,
                                       coordinates = reference$inference_object$table$index)
runtime <- proc.time()[[3L]] - started
saveRDS(precision, file.path(directory, "precision.rds"))
data <- prepare_analysis_data_v2(simulate_from_config_v2(request$config, request$seed),
                                  request$config, request$seed)$data
Sigma <- crossprod(data$X) / nrow(data$X)
rows <- reference$inference_object$table$index
one_step <- vapply(rows, function(k) {
  reference$fit_object$beta[k] - sum(precision$omega_symmetric[k, ] * reference$fit_object$score)
}, numeric(1))
variance <- vapply(rows, function(k) {
  omega <- precision$omega_symmetric[k, ]
  request$config$tau * (1 - request$config$tau) * drop(crossprod(omega, Sigma %*% omega)) / nrow(data$X)
}, numeric(1))
relative <- function(actual, expected) max(abs(actual - expected)) / max(1, max(abs(expected)))
old <- reference$inference_object$precision
checks <- data.frame(
  quantity = c("raw_precision", "symmetrised_precision", "one_step", "variance", "gamma", "all_original_LP_residuals"),
  discrepancy = c(relative(precision$omega_raw, old$omega_raw),
                    relative(precision$omega_symmetric, old$omega_symmetric),
                    relative(one_step, reference$inference_object$table$beta_tilde),
                    relative(variance, reference$inference_object$table$se^2),
                    abs(precision$gamma - old$gamma),
                    max(precision$matrix_grid_record$attempts$residual - precision$gamma)),
  tolerance = c(1e-5, 1e-5, 1e-5, 1e-5, 0, 1e-6)
)
checks$passed <- is.finite(checks$discrepancy) & checks$discrepancy <= checks$tolerance
utils::write.csv(checks, file.path(directory, "checks.csv"), row.names = FALSE)
metadata <- list(reference = normalizePath(reference_directory),
                   implementation_commit = current_commit_v2(root),
                   actual_LP_calculations = nrow(precision$matrix_grid_record$attempts),
                   reference_LP_calculations = nrow(old$matrix_grid_record$attempts),
                   runtime_sec = runtime,
                   certificate_max_residual = precision$matrix_grid_record$feasibility_certificate$max_residual,
                   minimum_gamma = precision$gamma,
                   certificate_verified = precision$matrix_grid_record$feasibility_certificate$verified,
                   nonselected_status = precision$raw_grid_records[[1L]]$status[-1L])
jsonlite::write_json(metadata, file.path(directory, "manifest.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
record <- precision$raw_grid_records[[1L]]$solver_records[[1L]]
writeLines(capture.output(str(record, max.level = 2L)), file.path(directory, "LP_record_fields.txt"))
print(checks)
if (!all(checks$passed) || !isTRUE(metadata$certificate_verified) || metadata$actual_LP_calculations != ncol(data$X)) {
  stop("The full-matrix certificate equivalence check failed.")
}
