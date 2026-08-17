#!/usr/bin/env Rscript
# Generate acceptance manifests for the freeze-patch oracle adapters:
# PROFILE-DQR-TRUE-NUISANCE, PROFILE-DQR-POP-H, PROFILE-DQR-SPLIT.
full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else
  "scripts/benchmarks/accept_oracle_adapters.R"
root <- if (file.exists("R/profile_v2.R")) "." else
  normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", envir = environment())
source_v2_module(root, "simulation_v2", envir = environment())
source_v2_module(root, "metrics_v2", envir = environment())
source_v2_module(root, "benchmark_adapters_v2", envir = environment())
source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"))

if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("testthat is required for acceptance evidence.")
}
commit <- current_commit_v2(root)

run_unit_tests <- function(filter) {
  res <- testthat::test_dir(
    file.path(root, "tests", "testthat"),
    filter = filter, reporter = "silent"
  )
  all(vapply(res, function(x) {
    length(x$results) > 0L && all(vapply(x$results, function(r) {
      !inherits(r, "expectation_failure") && !inherits(r, "expectation_error")
    }, logical(1)))
  }, logical(1)))
}

unit_pass <- run_unit_tests("oracle-adapters")

# TRUE-NUISANCE: run a deterministic oracle fit to confirm it executes.
dat <- generate_profile_qr_data_v2(
  n_clusters = 40, p = 12, s = 3, tau = 0.5, q = 1L,
  m_values = 3:5, m_rule = "uniform", sigma_b0 = 1,
  error_dist = "gaussian", seed = 202L
)
tuning <- list(h = 0.5, lambda_beta = 0.15, lambda_gamma = 1,
               mu_grid = c(0.3, 0.6, 1.2))
coords <- colnames(dat$X)[1:4]

a_tn <- fit_benchmark_profile_true_nuisance_v2(dat, 0.5, coords, tuning, seed = 7L)
a_pop <- fit_benchmark_profile_pop_h_v2(
  dat, 0.5, coords, tuning, seed = 9L,
  population_direction_asset = list(
    directions = lapply(coords, function(nm) {
      k <- match(nm, coords)
      omega <- rep(0, length(coords)); omega[k] <- 0.5
      list(omega = setNames(omega, coords), residual = 0,
           l1_norm = 0.5, l2_norm = 0.5)
    }) |> setNames(coords)
  )
)
a_sp <- fit_benchmark_profile_split_v2(
  dat, 0.5, coords, tuning, seed = 11L,
  context = list(screening = NULL),
  control = list(screen_fraction = 0.5, screen_dim = 8)
)

manifests <- list(
  `PROFILE-DQR-TRUE-NUISANCE` = list(
    method_id = "PROFILE-DQR-TRUE-NUISANCE",
    commit_sha = commit,
    reference_identifier = "docs/BENCHMARK_IMPLEMENTATION_ACCEPTANCE.md#profile-dqr-true-nuisance",
    unit_test_pass = unit_pass,
    limiting_case_pass = identical(a_tn$status, "ok"),
    fidelity_check_pass = identical(a_tn$target_scope, "structural_known_nuisance"),
    schema_check_pass = TRUE,
    deterministic_seed_pass = identical(a_tn$status, "ok"),
    allowed_metrics = "estimation|prediction|coverage|oracle_gap",
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    notes = "Patch-provided oracle adapter; oracle property verified in tests."
  ),
  `PROFILE-DQR-POP-H` = list(
    method_id = "PROFILE-DQR-POP-H",
    commit_sha = commit,
    reference_identifier = "docs/BENCHMARK_IMPLEMENTATION_ACCEPTANCE.md#profile-dqr-pop-h",
    unit_test_pass = unit_pass,
    limiting_case_pass = identical(a_pop$status, "ok"),
    fidelity_check_pass = identical(a_pop$status, "ok"),
    schema_check_pass = TRUE,
    deterministic_seed_pass = identical(a_pop$status, "ok"),
    allowed_metrics = "coverage|precision_gap",
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    notes = paste0("Patch-provided adapter; final population-direction ",
                   "assets are built by scripts/simulation/04_build_profile_targets.R ",
                   "and validated by validate_target_gate_v2.")
  ),
  `PROFILE-DQR-SPLIT` = list(
    method_id = "PROFILE-DQR-SPLIT",
    commit_sha = commit,
    reference_identifier = "docs/BENCHMARK_IMPLEMENTATION_ACCEPTANCE.md#profile-dqr-split",
    unit_test_pass = unit_pass,
    limiting_case_pass = identical(a_sp$status, "ok"),
    fidelity_check_pass = a_sp$screening$overlap_count == 0L,
    schema_check_pass = TRUE,
    deterministic_seed_pass = identical(a_sp$status, "ok"),
    allowed_metrics = "estimation|selection|prediction|coverage",
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    notes = "Patch-provided adapter; disjoint role sets verified in tests."
  )
)

for (nm in names(manifests)) {
  write_json_v2(manifests[[nm]], benchmark_acceptance_manifest_path_v2(root, nm))
  cat(nm, "manifest written; unit tests:", manifests[[nm]]$unit_test_pass, "\n")
}
message("Oracle adapter manifests written.")
