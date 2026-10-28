.libPaths(c(normalizePath("results/scratch/profile_pair_package_20261022/library"), .libPaths()))
library(jpDMEQR)
source("tests/testthat/helper-v2-skip.R")
Sys.setenv(JPDMEQR_CLUSTER_SOLVER = "compiled", JPDMEQR_ASSEMBLY = "batch",
           JPDMEQR_MATRIX_BACKEND = "Eigen_single_thread", JPDMEQR_IDENTITY_DETAIL = "aggregate",
           JPDMEQR_HESSIAN_REFRESH = "20")
result <- testthat::test_file("tests/testthat/test-profile-pair-dual-backend.R",
                              reporter = "summary", stop_on_failure = TRUE)
out <- "results/scratch/profile_dual_backend_20261022"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
saveRDS(result, file.path(out, "test_results.rds"))
summary <- as.data.frame(result)
write.csv(summary[, c("file", "test", "nb", "failed", "error", "skipped", "warning", "passed")],
          file.path(out, "test_summary.csv"), row.names = FALSE)
