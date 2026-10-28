root <- normalizePath(".", mustWork = TRUE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop("Supply certificate, grid, source, zero-start, strict-reference and new output run IDs.")
directories <- file.path(root, "results", "scratch", args)
output <- directories[6L]
if (dir.exists(output)) stop("The initial acceptance directory already exists.")
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", environment())
source_v2_module(root, "application_v2", environment())
actual <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
source_paths <- sort(c(
  file.path("R", list.files(file.path(root, "R"), "^v2_profile_v2_part[0-9]+[.]R$")),
  file.path("src", list.files(file.path(root, "src"), "[.](cpp|h|hpp)$", recursive = TRUE)),
  "src/Makevars", "src/Makevars.win", "DESCRIPTION", "NAMESPACE",
  file.path("R", sprintf("v2_application_v2_part%02d.R", c(3L, 5L, 6L, 9L, 10L, 13L))),
  "scripts/00_source_v2.R", "scripts/application/12_validate_GSE65391.R",
  "scripts/research/check_application_initial_certificate_20261023.R",
  "scripts/research/check_application_initial_clinical_20261023.R",
  "tests/testthat/test-application-initial-clinical-certificate.R",
  "config/application_execution_20261022.yml", "config/module_scoring_20261022.yml",
  "config/module_signatures.csv", "config/module_score_members_20261022.csv"
))
tracked <- system2("git", c("ls-files", "--error-unmatch", shQuote(source_paths)), stdout = FALSE, stderr = FALSE)
dirty <- system2("git", c("diff", "--quiet", "--ignore-cr-at-eol", "HEAD", "--", shQuote(source_paths)))
if (tracked != 0L || dirty != 0L) stop("The accepted initial numerical source must be committed and unchanged.")
read_checks <- function(path) {
  result <- read.csv(path, stringsAsFactors = FALSE)
  required <- c("nb", "failed", "error", "skipped", "warning", "passed")
  if (!nrow(result) || !all(required %in% names(result)) || anyNA(result[required]) ||
      any(result$failed != 0 | result$error | result$skipped | result$warning != 0) ||
      sum(result$nb) != sum(result$passed)) stop("A required actual initial numerical check failed.")
  list(path = normalizePath(path), assertions = sum(result$passed), pass = TRUE)
}
checks <- list(
  cache_and_full_design = read_checks(file.path(directories[1L], "unit_test_summary.csv")),
  fixed_penalty_stationarity_and_gene_continuation = read_checks(file.path(directories[1L], "test_summary.csv")),
  frozen_grid_global_kkt_and_timing = read_checks(file.path(directories[2L], "test_summary.csv")),
  matched_source_scope = read_checks(file.path(directories[3L], "test_summary.csv"))
)
manifest <- jsonlite::read_json(file.path(directories[1L], "completed_manifest.json"), simplifyVector = TRUE)
grid_manifest <- jsonlite::read_json(file.path(directories[2L], "completed_manifest.json"), simplifyVector = TRUE)
for (record in list(manifest, grid_manifest)) {
  ancestor <- system2("git", c("merge-base", "--is-ancestor", record$actual_source_commit, actual))
  changed <- system2("git", c("diff", "--quiet", "--ignore-cr-at-eol", record$actual_source_commit,
                               "--", shQuote(source_paths)))
  if (!identical(record$status, "completed") || ancestor != 0L || changed != 0L) {
    stop("The actual initial acceptance results must evaluate unchanged committed numerical source.")
  }
}
if (!identical(manifest$solver_configuration, grid_manifest$solver_configuration)) {
  stop("The initial certificate and timing grid must use the same actual solver configuration.")
}
grid <- read.csv(file.path(directories[2L], "scheduled_candidates.csv"), stringsAsFactors = FALSE)
timing <- read.csv(file.path(directories[2L], "actual_fit_summary.csv"), stringsAsFactors = FALSE)
differences <- read.csv(file.path(directories[2L], "actual_comparisons.csv"), stringsAsFactors = FALSE)
if (nrow(grid) != 108L || nrow(timing) != 216L || nrow(differences) != 108L ||
    !all(table(timing$candidate_id) == 2L) || !setequal(timing$mode, c("none", "clinical_kkt")) ||
    anyNA(timing[c("runtime_sec", "kkt_normalized")]) ||
    any(timing$runtime_sec < 0 | timing$kkt_normalized > 1e-3) ||
    any(!timing$final_certificate[timing$mode == "clinical_kkt"]) ||
    any(!differences$selection_identical)) stop("The complete actual 108-candidate numerical grid is required.")
fixed <- read.csv(file.path(directories[1L], "fixed_penalty_stationarity_summary.csv"), stringsAsFactors = FALSE)
if (nrow(fixed) != 6L || !setequal(fixed$q, 1:2) || !setequal(fixed$tau, c(.25, .5, .75))) {
  stop("All six actual fixed-penalty stationarity checks are required.")
}
low <- readRDS(file.path(directories[1L], "actual_low_penalty_continuation.rds"))
low_selected <- sum(abs(low$beta_l1[low$penalty_factor > 0]) > 1e-10)
if (low_selected < 1L || !isTRUE(low$converged) || low$kkt_normalized > 1e-3 ||
    any(vapply(low$initial_optimizer_history, function(pass) pass$global_certificate_accepted, logical(1)))) {
  stop("The actual nonzero-gene continuation branch must be observed and accepted.")
}
zero <- read.csv(file.path(directories[4L], "actual_comparisons.csv"), stringsAsFactors = FALSE)
strict <- readRDS(file.path(directories[5L], "failure_condition.rds"))
if (nrow(zero) != 6L || !inherits(strict$condition, "condition")) {
  stop("The saved zero-start stopping diagnostics must be available.")
}
paths <- list(
  certificate_manifest = file.path(directories[1L], "completed_manifest.json"),
  fixed_penalty_stationarity = file.path(directories[1L], "fixed_penalty_stationarity_summary.csv"),
  nonzero_gene_continuation = file.path(directories[1L], "actual_low_penalty_continuation.rds"),
  grid_manifest = file.path(directories[2L], "completed_manifest.json"),
  grid_timing = file.path(directories[2L], "actual_fit_summary.csv"),
  grid_differences = file.path(directories[2L], "actual_comparisons.csv"),
  zero_start_stopping_differences = file.path(directories[4L], "actual_comparisons.csv"),
  strict_zero_start_condition = file.path(directories[5L], "failure_condition.rds"),
  matched_call_graph = file.path(directories[3L], "actual_matched_call_graph.rds"),
  matched_source_identity = file.path(directories[3L], "actual_source_identity.rds")
)
if (any(!file.exists(unlist(paths)))) stop("Every actual acceptance source path is required.")
dir.create(output, recursive = TRUE)
result <- list(
  pass = TRUE, method_version = "profile_pair_20261022_initial_l1",
  actual_source_commit = actual, source_paths = source_paths,
  acceptance_builder = list(path = "scripts/research/freeze_application_initial_acceptance_20261023.R",
                             actual_source_commit = actual,
                             scope = "metadata_validation_without_beta_or_profile_score_calculation"),
  initial_configuration = list(clinical_acceleration = "clinical_kkt"),
  solver_configuration = manifest$solver_configuration,
  checks = checks, actual_check_paths = paths,
  fixed_penalty_scope = "stationarity_at_globally_certified_point",
  zero_start_diagnostics = list(
    scope = "development_reference_stopping_diagnostics",
    strict_coefficient_equivalence_pass = FALSE, required_for_global_kkt_certificate = FALSE,
    max_coefficient_difference = max(zero$beta_difference),
    max_prediction_difference = max(zero$prediction_difference),
    max_objective_difference = max(zero$objective_difference),
    strict_reference_failure_message = conditionMessage(strict$condition),
    strict_reference_kkt_normalized = strict$condition$diagnostics$kkt_normalized
  ),
  nonzero_gene_continuation = list(pass = TRUE, selected_genes = low_selected,
                                   kkt_normalized = low$kkt_normalized, research_lambda_multiplier = .01),
  frozen_grid = list(pass = TRUE, candidates = 108L, actual_fits = 216L,
                      original_runtime_sec = sum(timing$runtime_sec[timing$mode == "none"]),
                      accelerated_runtime_sec = sum(timing$runtime_sec[timing$mode == "clinical_kkt"]),
                      strict_zero_start_differences_recorded = TRUE,
                      all_selection_indicators_identical = all(differences$selection_identical)),
  completed_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
)
jsonlite::write_json(result, file.path(output, "initial_acceptance_manifest.json"),
                     auto_unbox = TRUE, pretty = TRUE, digits = NA)
write.csv(data.frame(check = names(checks), assertions = vapply(checks, `[[`, numeric(1), "assertions"),
                       pass = TRUE), file.path(output, "acceptance_check_counts.csv"), row.names = FALSE)
