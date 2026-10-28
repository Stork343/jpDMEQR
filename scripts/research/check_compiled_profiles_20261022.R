.libPaths(c(normalizePath("results/scratch/profile_pair_package_20261022/library"), .libPaths()))
library(jpDMEQR)
source("tests/testthat/helper-v2-skip.R")
result <- testthat::test_file("tests/testthat/test-compiled-profile-parity.R",
                              reporter = "summary", stop_on_failure = TRUE)
out <- "results/scratch/compiled_profiles_20261022"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
saveRDS(result, file.path(out, "test_results.rds"))
summary <- as.data.frame(result)
write.csv(summary[, c("file", "test", "nb", "failed", "error", "skipped", "warning", "passed")],
          file.path(out, "test_summary.csv"), row.names = FALSE)
stopifnot(!any(summary$failed), !any(summary$error), !any(summary$skipped))
