args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1]] else ".", mustWork = TRUE)
setwd(root)
library_path <- file.path(root, "results", "scratch", "profile_pair_package_20261022", "library")
if (dir.exists(library_path)) .libPaths(c(normalizePath(library_path), .libPaths()))
output <- if (length(args) >= 2L) args[[2]] else
  "results/scratch/profile_pair_tuning_20261022/actual_dgp_01"
if (dir.exists(output)) stop("The tuning audit output directory already exists.")
dir.create(output, recursive = TRUE)
source(file.path("tests", "testthat", "helper-v2-skip.R"))
test_environment <- new.env(parent = globalenv())
tuning_records <- new.env(parent = emptyenv())
options(jpDMEQR.frozen_tuning_test_results = tuning_records)
results <- testthat::test_file(
  file.path("tests", "testthat", "test-frozen-penalty-tuning.R"),
  reporter = "summary", env = test_environment, stop_on_failure = FALSE
)
table <- as.data.frame(results)
saveRDS(results, file.path(output, "test_results.rds"))
utils::write.csv(table[, c("file", "test", "failed", "error", "warning", "skipped", "passed")],
                 file.path(output, "assertions.csv"), row.names = FALSE)
records <- as.list(tuning_records)
saveRDS(records, file.path(output, "actual_tuning_records.rds"))
required <- c("paired", "oracle_support", "POOLED-QR-LASSO", "SQR-DEBIASED-IID",
              "POOLED-QR-LASSO_full", "SQR-DEBIASED-IID_full", "failed_candidates", "runner")
record_complete <- all(required %in% names(records))
pass <- record_complete && !any(table$failed > 0 | table$error | table$warning > 0 | table$skipped)
jsonlite::write_json(list(
  pass = pass, assertions = sum(table$passed), failures = sum(table$failed),
  errors = sum(table$error), warnings = sum(table$warning), skipped = sum(table$skipped),
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  data_source = "generate_profile_qr_data_v2_from_registered_E07_design",
  method_version = "profile_pair_20261022",
  records_complete = record_complete, record_names = names(records)
), file.path(output, "manifest.json"), auto_unbox = TRUE, pretty = TRUE)
writeLines(capture.output(sessionInfo()), file.path(output, "sessionInfo.txt"))
if (!pass) stop("The frozen tuning audit failed; inspect ", output)
writeLines("complete", file.path(output, "COMPLETE"))
