args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Expected repository, test file and output directory.")
root <- normalizePath(args[[1]], mustWork = TRUE)
setwd(root)
test_file <- normalizePath(args[[2]], mustWork = TRUE)
output <- normalizePath(args[[3]], mustWork = TRUE)
library_path <- normalizePath("results/scratch/profile_pair_package_20261022/library", mustWork = TRUE)
.libPaths(c(library_path, .libPaths()))
library(jpDMEQR)
library(testthat)
if (!identical(as.character(packageVersion("jpDMEQR")), "0.4.0") ||
    !is.loaded("jpdmeqr_cluster_profiles", PACKAGE = "jpDMEQR") ||
    !is.loaded("jpDMEQR_population_crossprod", PACKAGE = "jpDMEQR")) {
  stop("The installed 0.4.0 package and both registered native functions are required.")
}
configuration <- c(
  JPDMEQR_CLUSTER_SOLVER = "compiled", JPDMEQR_ASSEMBLY = "batch",
  JPDMEQR_MATRIX_BACKEND = "Eigen_single_thread", JPDMEQR_IDENTITY_DETAIL = "aggregate",
  JPDMEQR_UNPENALISED_OPTIMIZER = "newton_armijo", JPDMEQR_PROXIMAL_SCALING = "cluster_rms",
  JPDMEQR_DANTZIG_SOLVER = "dual_warm", JPDMEQR_HESSIAN_REFRESH = "20",
  JPDMEQR_SHARED_FIRST_STAGE = "TRUE", JPDMEQR_CV_CORES = "1",
  JPDMEQR_INITIAL_ACCELERATION = if (basename(test_file) == "test-application-initial-estimator.R")
    "none" else "clinical_kkt", NOT_CRAN = "true"
)
do.call(Sys.setenv, as.list(configuration))
source(file.path("tests", "testthat", "helper-v2-skip.R"))
testthat::set_max_fails(Inf)
options(jpDMEQR.unit_artifact_dir = output)
started <- Sys.time()
source_error <- NULL
results <- tryCatch(
  testthat::test_file(test_file, reporter = "summary", stop_on_failure = FALSE),
  error = function(condition) {
    source_error <<- condition
    NULL
  }
)
saveRDS(list(results = results, source_error = source_error), file.path(output, "test_results.rds"))
table <- if (is.null(results)) data.frame() else as.data.frame(results)
if (nrow(table)) {
  utils::write.csv(table[, c("file", "test", "nb", "failed", "error", "warning", "skipped", "passed")],
                   file.path(output, "assertions.csv"), row.names = FALSE)
}
counts <- list(passed = sum(table$passed), failed = sum(table$failed),
               errors = sum(table$error), warnings = sum(table$warning),
               skipped = sum(table$skipped), test_blocks = nrow(table))
warning_rows <- list()
for (block in results) {
  for (condition in block$results) {
    if (inherits(condition, "expectation_warning")) {
      warning_rows[[length(warning_rows) + 1L]] <- data.frame(
        file = basename(test_file), test = block$test,
        line = if (is.null(condition$srcref)) NA_integer_ else as.integer(condition$srcref[[1L]]),
        message = conditionMessage(condition), stringsAsFactors = FALSE
      )
    }
  }
}
warning_details <- if (length(warning_rows)) do.call(rbind, warning_rows) else
  data.frame(file = character(), test = character(), line = integer(), message = character())
utils::write.csv(warning_details, file.path(output, "warning_details.csv"), row.names = FALSE)
if (nrow(warning_details) != counts$warnings) stop("Warning records do not reproduce the test count.")
assertion_pass <- is.null(source_error) && nrow(table) > 0L &&
  counts$failed == 0L && counts$errors == 0L && counts$skipped == 0L
strict_pass <- assertion_pass && counts$warnings == 0L
manifest <- c(list(
  pass = assertion_pass, assertion_pass = assertion_pass, strict_pass = strict_pass,
  warning_count = counts$warnings, file = basename(test_file), configuration = as.list(configuration),
  source_error = if (is.null(source_error)) NULL else conditionMessage(source_error),
  package_version = as.character(packageVersion("jpDMEQR")),
  installed_library = library_path, temp_directory = tempdir(),
  locale = Sys.getlocale(),
  started_utc = format(started, tz = "UTC", usetz = TRUE),
  elapsed_sec = as.numeric(difftime(Sys.time(), started, units = "secs"))
), counts)
jsonlite::write_json(manifest, file.path(output, "manifest.json"),
                     auto_unbox = TRUE, pretty = TRUE, null = "null")
writeLines(capture.output(sessionInfo()), file.path(output, "sessionInfo.txt"))
if (!assertion_pass) quit(status = 1L)
