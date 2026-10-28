root <- normalizePath(".", mustWork = TRUE)
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2")
source_v2_module(root, "simulation_v2")
source(file.path(root, "tests", "testthat", "helper-v2-skip.R"))
out <- file.path(root, "results", "scratch", "profile_pair_population_performance_20261022")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(JPDMEQR_CLUSTER_SOLVER = "newton")
mode <- Sys.getenv("JPDMEQR_POPULATION_TEST_BACKEND", unset = "R_BLAS")
if (mode %in% c("Eigen", "compiled_Eigen")) {
  .libPaths(c(normalizePath(file.path(root, "results", "scratch", "profile_pair_package_20261022", "library")), .libPaths()))
  library(jpDMEQR)
  profile_pair_population_eigen_init_v2(root)
  if (mode == "compiled_Eigen") options(jpDMEQR.population_nuisance_backend = "compiled")
}
out <- file.path(out, mode)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
result <- testthat::test_dir(file.path(root, "tests", "testthat"),
                              filter = "profile-pair-((streamed|conditional)-population|population-eigen)",
                              reporter = "summary", stop_on_failure = FALSE)
saveRDS(result, file.path(out, "equivalence_tests.rds"))
summary <- as.data.frame(result)
utils::write.csv(summary[, c("file", "test", "nb", "failed", "error", "skipped", "warning", "passed")],
                 file.path(out, "equivalence_tests.csv"), row.names = FALSE)
stopifnot(!any(summary$error), all(summary$failed == 0L), sum(summary$skipped) == 0L)
print(c(assertions = sum(summary$nb), passed = sum(summary$passed)))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
