args <- base::commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1]] else ".", mustWork = TRUE)
setwd(root)
source(file.path("tests", "testthat", "helper-v2-skip.R"))
result <- testthat::test_dir(
  file.path("tests", "testthat"),
  filter = "profile-(score-repeat-limit|bandwidth-pair)",
  reporter = "summary", stop_on_failure = FALSE
)
saveRDS(result, file.path("results", "scratch", "profile_score_audit_20261022",
                         "testthat_results.rds"))
summary <- as.data.frame(result)
utils::write.csv(summary[, c("file", "test", "nb", "failed", "error", "skipped", "warning", "passed")],
                 file.path("results", "scratch", "profile_score_audit_20261022",
                           "testthat_summary.csv"), row.names = FALSE)
stopifnot(all(summary$failed == 0), !any(summary$error),
          sum(summary$skipped) == 0L)
