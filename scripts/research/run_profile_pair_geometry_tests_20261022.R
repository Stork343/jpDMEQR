root <- normalizePath(".", mustWork = TRUE)
source(file.path(root, "tests", "testthat", "helper-v2-skip.R"))
out <- file.path(root, "results", "scratch", "profile_pair_geometry_tests_20261022")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
result <- testthat::test_file(
  file.path(root, "tests", "testthat", "test-profile-pair-geometry-audit.R"),
  reporter = "summary", stop_on_failure = FALSE
)
saveRDS(result, file.path(out, "results.rds"))
summary <- as.data.frame(result)
utils::write.csv(summary[, c("file", "test", "nb", "failed", "error", "skipped", "warning", "passed")],
                 file.path(out, "summary.csv"), row.names = FALSE)
stopifnot(!any(summary$error), all(summary$failed == 0L))
print(c(assertions = sum(summary$nb), passed = sum(summary$passed),
        failed = sum(summary$failed), skipped = sum(summary$skipped)))
