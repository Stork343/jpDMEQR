args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[1] else ".", mustWork = TRUE)
.libPaths(c(file.path(root, "results", "scratch", "profile_pair_package_20261022", "library"), .libPaths()))
library(jpDMEQR)
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", environment())
source(file.path(root, "tests", "testthat", "helper-v2-skip.R"))
configuration <- solver_configuration_v2()
if (configuration$cluster_solver != "compiled" || configuration$proximal_scaling != "cluster_rms" ||
    configuration$unpenalised_optimizer != "newton_armijo") stop("The complete production numerical configuration is required.")
directory <- file.path(root, "results", "scratch", "application_exploratory_tests_20261022",
                        paste0("production_", format(Sys.time(), "%Y%m%dT%H%M%OS6")))
dir.create(directory, recursive = TRUE)
result <- testthat::test_file(file.path(root, "tests", "testthat", "test-application-exploratory-pair.R"),
                               reporter = "summary", stop_on_failure = FALSE)
saveRDS(result, file.path(directory, "test_results.rds"))
table <- as.data.frame(result)
table$result <- NULL
utils::write.csv(table, file.path(directory, "test_summary.csv"), row.names = FALSE)
passed <- nrow(table) > 0L && all(table$failed == 0L) && !any(table$error) && !any(table$skipped)
jsonlite::write_json(list(commit_sha = system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"), stdout = TRUE),
                          solver_configuration = configuration, unit_checks_pass = passed,
                          assertions = sum(table$passed), warnings = sum(table$warning),
                          source_test = "tests/testthat/test-application-exploratory-pair.R",
                          method_version = "profile_pair_20261022"),
                     file.path(directory, "manifest.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
if (!passed) stop("The production exploratory source checks failed.")
message("Exploratory source checks saved to ", directory)
