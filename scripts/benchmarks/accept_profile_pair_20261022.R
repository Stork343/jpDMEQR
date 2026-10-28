args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[1] else ".", mustWork = TRUE)
.libPaths(c(file.path(root, "results", "scratch", "profile_pair_package_20261022", "library"), .libPaths()))
library(jpDMEQR)
destination <- if (length(args) >= 2L) args[2] else
  file.path(root, "results", "preflight", "benchmarks_profile_pair_20261022")
source(file.path(root, "scripts", "00_source_v2.R"))
for (module in c("profile_v2", "simulation_v2", "metrics_v2", "benchmark_adapters_v2")) {
  source_v2_module(root, module, environment())
}
source(file.path(root, "tests", "testthat", "helper-v2-skip.R"))
commit <- system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"), stdout = TRUE)
if (length(commit) != 1L || !grepl("^[0-9a-f]{40}$", commit)) {
  stop("Acceptance requires the actual repository commit identity.")
}
files <- c("test-benchmark-profile-pair.R", "test-profile-bandwidth-pair.R",
           "test-profile-pair-fit-inference.R", "test-profile-pair-population-assets.R",
           "test-profile-pair-dual-backend.R")
results <- lapply(files, function(path) {
  testthat::test_file(file.path(root, "tests", "testthat", path),
                      reporter = "summary", stop_on_failure = FALSE)
})
names(results) <- files
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
saveRDS(results, file.path(destination, "test_results.rds"))
test_pass <- function(test) {
  length(test$results) > 0L && !any(vapply(test$results, function(result) {
    inherits(result, c("expectation_failure", "expectation_error", "expectation_skip"))
  }, logical(1)))
}
core_pass <- all(vapply(results[-1L], function(result) {
  length(result) > 0L && all(vapply(result, test_pass, logical(1)))
}, logical(1)))
solver_configuration <- solver_configuration_v2()
solver_configuration$paired_dantzig_backend <- profile_pair_precision_backend_v2()
solver_configuration$numerical_reference <- list(backend = "clarabel_native",
                                                  static_regularization_constant = 1e-6)
solver_configuration$dual_reference <- list(backend = "dual_warm", solver = "HiGHS-dual-simplex",
                                            threads = 1L, presolve = "off",
                                            primal_feasibility_tolerance = 1e-8,
                                            dual_feasibility_tolerance = 1e-8)
benchmark_tests <- results[[1L]]
methods <- c("PROFILE-DQR", "PROFILE-DQR-TRUE-SUPPORT", "PROFILE-DQR-TRUE-NUISANCE",
             "PROFILE-DQR-POP-H", "PROFILE-DQR-SPLIT", "POSTREFIT-EXACT-H")
metric_permissions <- c(
  "estimation|selection|prediction|coverage|theory|computation",
  "estimation|prediction|coverage|oracle_gap",
  "estimation|prediction|coverage|oracle_gap", "coverage|precision_gap",
  "estimation|selection|prediction|coverage", "estimation|precision_gap")
rows <- vector("list", length(methods))
for (i in seq_along(methods)) {
  method <- methods[i]
  relevant <- benchmark_tests[vapply(benchmark_tests, function(test) {
    startsWith(test$test, paste(method, "paired"))
  }, logical(1))]
  check_category <- function(category) {
    subset <- relevant[vapply(relevant, function(test) {
      grepl(category, test$test, fixed = TRUE)
    }, logical(1))]
    length(subset) > 0L && all(vapply(subset, test_pass, logical(1)))
  }
  manifest <- list(
    method_id = method, method_version = "profile_pair_20261022",
    score_definition = "profile_bandwidth_pair_2_half_minus_full",
    solver_configuration = solver_configuration,
    commit_sha = commit, reference_identifier = "docs/METHOD_SPECIFICATION.md",
    unit_test_pass = length(relevant) > 0L && all(vapply(relevant, test_pass, logical(1))),
    limiting_case_pass = check_category("limiting case"),
    fidelity_check_pass = core_pass && check_category("fidelity"),
    schema_check_pass = check_category("schema"),
    deterministic_seed_pass = check_category("deterministic"),
    allowed_metrics = metric_permissions[i],
    test_names = vapply(relevant, `[[`, character(1), "test"),
    assertion_count = sum(vapply(relevant, function(test) sum(vapply(test$results, function(value) {
      inherits(value, "expectation_success")
    }, logical(1))), integer(1))),
    warning_count = sum(vapply(relevant, function(test) sum(vapply(test$results, function(value) {
      inherits(value, "expectation_warning")
    }, logical(1))), integer(1))),
    population_fixture = if (method == "PROFILE-DQR-POP-H")
      "independent genuine DGP development assets from the paired population builders" else NULL,
    final_population_assets_accepted = FALSE,
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
  flags <- unlist(manifest[c("unit_test_pass", "limiting_case_pass", "fidelity_check_pass",
                            "schema_check_pass", "deterministic_seed_pass")])
  path <- file.path(destination, method)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(manifest, file.path(path, "acceptance.json"),
                        auto_unbox = TRUE, pretty = TRUE, null = "null")
  rows[[i]] <- data.frame(method_id = method, passed = all(flags),
                           assertions = manifest$assertion_count,
                           stringsAsFactors = FALSE)
}
summary <- do.call(rbind, rows)
utils::write.csv(summary, file.path(destination, "acceptance_summary.csv"), row.names = FALSE)
print(summary)
if (!all(summary$passed) || !core_pass) stop("Paired benchmark acceptance failed.")
