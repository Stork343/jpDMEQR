source(file.path("scripts", "00_source_v2.R"))
source(file.path("tests", "testthat", "helper-v2-skip.R"))
for (module in c("profile_v2", "simulation_v2", "metrics_v2", "benchmark_adapters_v2",
                 "application_v2")) {
  source_v2_module(".", module, .GlobalEnv)
}
Sys.setenv(JPDMEQR_CLUSTER_SOLVER = "newton", JPDMEQR_DANTZIG_SOLVER = "clarabel_native",
           JPDMEQR_ASSEMBLY = "sequential", JPDMEQR_CV_CORES = "1",
           JPDMEQR_HESSIAN_REFRESH = "1", JPDMEQR_SHARED_FIRST_STAGE = "FALSE")
results <- list()
for (path in c("test-profile-bandwidth-pair.R", "test-profile-pair-fit-inference.R",
               "test-profile-pair-population-assets.R", "test-profile-pair-geometry-audit.R",
               "test-profile-pair-runner-accounting.R", "test-application-module-sources.R",
               "test-application-profile-pair.R", "test-application-profile-workflow.R")) {
  results[[path]] <- testthat::test_file(
    file.path("tests", "testthat", path), reporter = "summary", stop_on_failure = TRUE
  )
}
snapshot <- memory_snapshot_v2()
stopifnot(is.finite(snapshot[["rss_mb"]]), snapshot[["rss_mb"]] > 0)
out_dir <- file.path("results", "scratch", "profile_pair_integration_20261022")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(results, file.path(out_dir, "test_results.rds"))
summary <- do.call(rbind, lapply(results, as.data.frame))
stopifnot(all(summary$failed == 0), !any(summary$error), !any(summary$skipped))
write.csv(summary[, c("file", "test", "nb", "failed", "error", "skipped", "warning", "passed")],
          file.path(out_dir, "test_summary.csv"), row.names = FALSE)
write.csv(data.frame(metric = names(snapshot), value = unname(snapshot)),
          file.path(out_dir, "memory_snapshot.csv"), row.names = FALSE)
message("Paired method integration checks passed.")
