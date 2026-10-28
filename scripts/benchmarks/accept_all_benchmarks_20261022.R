run_benchmark_acceptance_20261022 <- function(root = ".", methods = NULL,
                                               reference_directory = NULL) {
  root <- normalizePath(root, mustWork = TRUE)
  .libPaths(c(file.path(root, "results", "scratch", "profile_pair_package_20261022", "library"), .libPaths()))
  library(jpDMEQR)
  source(file.path(root, "scripts", "00_source_v2.R"), local = TRUE)
  for (module in c("profile_v2", "simulation_v2", "metrics_v2", "benchmark_adapters_v2")) {
    source_v2_module(root, module, environment())
  }
  source(file.path(root, "tests", "testthat", "helper-v2-skip.R"), local = FALSE)
  requirements <- utils::read.csv(file.path(root, "config", "benchmark_requirements.csv"),
                                  stringsAsFactors = FALSE)
  baseline <- c("POOLED-QR-LASSO", "LQMM", "SQR-DEBIASED-IID", "BIAS-ADJ-LQMM",
                "DOUBLE-PEN-QLMM", "QGEE-SCAD", "QIF-SEE")
  paired <- c("PROFILE-DQR", "PROFILE-DQR-TRUE-SUPPORT", "PROFILE-DQR-TRUE-NUISANCE",
              "PROFILE-DQR-POP-H", "PROFILE-DQR-SPLIT", "POSTREFIT-EXACT-H")
  if (is.null(methods)) methods <- c(paired, baseline)
  if (length(setdiff(methods, c(paired, baseline)))) stop("Unknown acceptance method.")
  commit <- system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"), stdout = TRUE)
  if (length(commit) != 1L) stop("The actual implementation commit is required.")
  destination <- file.path(root, "results", "preflight", "benchmarks")
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  timestamp <- format(Sys.time(), "%Y%m%dT%H%M%OS6", tz = "UTC")
  records <- file.path(root, "results", "preflight", paste0("benchmark_acceptance_", timestamp))
  if (dir.exists(records)) stop("The acceptance record directory already exists.")
  dir.create(records, recursive = TRUE)
  rscript <- file.path(R.home("bin"), "Rscript")
  if (length(intersect(methods, paired))) {
    status <- system2(rscript, c(shQuote(file.path(root, "scripts", "benchmarks", "accept_profile_pair_20261022.R")),
                                  shQuote(root), shQuote(destination)))
    if (status != 0L) stop("The actual paired acceptance suite failed.")
  }
  selected_baseline <- intersect(methods, baseline)
  if (!length(selected_baseline)) return(invisible(records))
  if (is.null(reference_directory)) {
    reference_directory <- file.path(records, "reference_fidelity")
    status <- system2(rscript, c(shQuote(file.path(root, "scripts", "benchmarks", "check_reference_fidelity_20261022.R")),
                                  shQuote(root), shQuote(reference_directory)))
    if (status != 0L) stop("The actual source-fidelity calculation failed.")
  }
  reference <- utils::read.csv(file.path(reference_directory, "checks.csv"), stringsAsFactors = FALSE)
  reference_manifest <- jsonlite::read_json(file.path(reference_directory, "reference_fidelity.json"), simplifyVector = TRUE)
  if (!identical(reference_manifest$commit_sha, commit)) stop("Source-fidelity checks require the current commit.")
  files <- unique(c("test-benchmark-patch-adapters.R", "test-benchmark-quantile-l1-block.R",
                    paste0("test-benchmark-", setdiff(selected_baseline, c("POOLED-QR-LASSO", "LQMM")), ".R")))
  if ("SQR-DEBIASED-IID" %in% selected_baseline) files <- c(files, "test-benchmark-SQR-CV.R")
  results <- lapply(files, function(file) {
    testthat::test_file(file.path(root, "tests", "testthat", file), reporter = "summary", stop_on_failure = FALSE)
  })
  names(results) <- files
  saveRDS(results, file.path(records, "test_results.rds"))
  tests <- unlist(results, recursive = FALSE)
  test_pass <- function(test) {
    length(test$results) > 0L && !any(vapply(test$results, function(value) {
      inherits(value, c("expectation_failure", "expectation_error", "expectation_skip"))
    }, logical(1)))
  }
  category <- function(relevant, pattern) {
    matches <- relevant[vapply(relevant, function(test) grepl(pattern, test$test), logical(1))]
    length(matches) > 0L && all(vapply(matches, test_pass, logical(1)))
  }
  rows <- list()
  for (method in selected_baseline) {
    relevant <- tests[vapply(tests, function(test) startsWith(test$test, method), logical(1))]
    if (method == "POOLED-QR-LASSO") {
      block <- results[["test-benchmark-quantile-l1-block.R"]]
      relevant <- c(relevant, block)
    }
    source_checks <- reference[reference$method_id == method, , drop = FALSE]
    source_pass <- nrow(source_checks) > 0L && all(source_checks$passed) &&
      all(is.finite(source_checks$discrepancy) & source_checks$discrepancy <= source_checks$tolerance)
    limiting_pattern <- switch(method,
      "POOLED-QR-LASSO" = "scalar random intercept blocks",
      "LQMM" = "limiting case", "SQR-DEBIASED-IID" = "limiting case",
      "BIAS-ADJ-LQMM" = "variance limit", "DOUBLE-PEN-QLMM" = "suppressed nuisance",
      "QGEE-SCAD" = "independence working correlation", "QIF-SEE" = "source limiting case")
    rr <- requirements[requirements$method_id == method, , drop = FALSE]
    manifest <- list(
      method_id = method, commit_sha = commit, reference_identifier = rr$reference_identifier,
      method_version = if (method == "SQR-DEBIASED-IID") "source_20261023_CV" else NULL,
      unit_test_pass = length(relevant) > 0L && all(vapply(relevant, test_pass, logical(1))),
      limiting_case_pass = category(relevant, limiting_pattern),
      fidelity_check_pass = source_pass,
      schema_check_pass = category(relevant, "schema|deterministic and pooled|low-dimensional set"),
      deterministic_seed_pass = category(relevant, "deterministic"),
      allowed_metrics = rr$allowed_metrics,
      source_checks = source_checks, reference_record = normalizePath(reference_directory),
      test_record = normalizePath(file.path(records, "test_results.rds")),
      test_names = vapply(relevant, `[[`, character(1), "test"),
      assertion_count = sum(vapply(relevant, function(test) sum(vapply(test$results, function(value) {
        inherits(value, "expectation_success")
      }, logical(1))), integer(1))),
      warning_count = sum(vapply(relevant, function(test) sum(vapply(test$results, function(value) {
        inherits(value, "expectation_warning")
      }, logical(1))), integer(1))),
      source_specification_status = if (method == "QIF-SEE")
        "Algorithm1_explicit_density_and_ordinary_inverse_domain" else "complete_source_or_maintained_package_checked",
      created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
    )
    path <- file.path(destination, method)
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    jsonlite::write_json(manifest, file.path(path, "acceptance.json"),
                          auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
    flags <- unlist(manifest[c("unit_test_pass", "limiting_case_pass", "fidelity_check_pass",
                              "schema_check_pass", "deterministic_seed_pass")])
    rows[[length(rows) + 1L]] <- data.frame(method_id = method, passed = all(flags),
                                            assertions = manifest$assertion_count)
  }
  table <- do.call(rbind, rows)
  utils::write.csv(table, file.path(records, "acceptance_summary.csv"), row.names = FALSE)
  print(table)
  if (!all(table$passed)) stop("One or more benchmark acceptance requirements remain incomplete.")
  invisible(records)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  run_benchmark_acceptance_20261022(if (length(args)) args[1] else ".")
}
