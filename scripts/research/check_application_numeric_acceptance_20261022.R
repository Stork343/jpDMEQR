root <- normalizePath(".", mustWork = TRUE)
.libPaths(c(file.path(root, "results", "scratch", "profile_pair_package_20261022", "library"), .libPaths()))
library(jpDMEQR)
Sys.setenv(JPDMEQR_CLUSTER_SOLVER = "compiled",
           JPDMEQR_UNPENALISED_OPTIMIZER = "newton_armijo",
           JPDMEQR_PROXIMAL_SCALING = "cluster_rms")
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", environment())
source_v2_module(root, "simulation_v2", environment())
source_v2_module(root, "application_v2", environment())
v2_require_root <- function() root
args <- commandArgs(trailingOnly = TRUE)
run_id <- if (length(args)) args[1] else "numerical_solver_acceptance_20261022_v3"
output <- file.path(root, "results", "scratch", run_id)
if (dir.exists(output)) stop("The numerical acceptance directory already exists.")
dir.create(output, recursive = TRUE)
files <- c("test-profile-unpenalised-newton.R", "test-profile-proximal-scaling.R",
           "test-application-initial-estimator.R", "test-application-compact-workflow.R")
results <- lapply(files, function(file) {
  as.data.frame(testthat::test_file(file.path(root, "tests", "testthat", file), reporter = "summary"))
})
table <- do.call(rbind, results)
fields <- c("file", "context", "test", "nb", "failed", "skipped", "error", "warning", "passed")
write.csv(table[fields], file.path(output, "test_summary.csv"), row.names = FALSE)
jsonlite::write_json(solver_configuration_v2(), file.path(output, "solver_configuration.json"),
                     auto_unbox = TRUE, pretty = TRUE)
manifest <- list(
  implementation_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  dependency_paths = application_numeric_dependency_paths_v2(root),
  assertions = sum(table$nb), passed_assertions = sum(table$passed),
  failed_assertions = sum(table$failed), error_cases = sum(table$error),
  skipped_cases = sum(table$skipped), warning_cases = sum(table$warning),
  solver_configuration = solver_configuration_v2(),
  completed_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
)
jsonlite::write_json(manifest, file.path(output, "validation_manifest.json"), auto_unbox = TRUE, pretty = TRUE)
if (sum(table$passed) != 148L || sum(table$nb) != 148L ||
    any(table$failed > 0 | table$error | table$skipped | table$warning > 0)) {
  stop("Production numerical acceptance failed.")
}
