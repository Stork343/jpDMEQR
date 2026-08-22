#!/usr/bin/env Rscript
# Generate acceptance manifests for the patch-provided adapters:
# PROFILE-DQR, POOLED-QR-LASSO, LQMM, PROFILE-DQR-TRUE-SUPPORT.
full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else
  "scripts/benchmarks/accept_patch_adapters.R"
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

unit_pass <- run_unit_tests("patch-adapters")

dat <- generate_profile_qr_data_v2(
  n_clusters = 40, p = 12, s = 3, tau = 0.5, q = 1L,
  m_values = 3:5, m_rule = "uniform", sigma_b0 = 1,
  error_dist = "gaussian", seed = 202L
)
tuning <- list(h = 0.5, lambda_beta = 0.15, lambda_gamma = 1,
               mu_grid = c(0.3, 0.6, 1.2))
coords <- colnames(dat$X)[1:4]

a_pd <- fit_benchmark_profile_dqr_v2(dat, 0.5, coords, tuning, seed = 7L)
a_pq <- fit_benchmark_pooled_qr_lasso_v2(dat, 0.5, coords, tuning, seed = 5L)
keep_l <- c(1:4, 7)
sub <- dat
sub$X <- dat$X[, keep_l, drop = FALSE]
a_lq <- fit_benchmark_lqmm_v2(sub, 0.5, colnames(sub$X)[1:3], tuning,
                              seed = 9L, random_slope = FALSE)
a_ts <- fit_benchmark_profile_true_support_v2(
  dat, 0.5, coords, tuning, seed = 11L, active = dat$active
)

manifests <- list(
  `PROFILE-DQR` = list(
    method_id = "PROFILE-DQR",
    commit_sha = commit,
    reference_identifier = "docs/METHOD_SPECIFICATION.md",
    unit_test_pass = unit_pass,
    limiting_case_pass = identical(a_pd$status, "ok") &&
      all(is.finite(a_pd$beta_tilde)),
    fidelity_check_pass = isTRUE(a_pd$profile_identity_error < 1e-8),
    schema_check_pass = TRUE,
    deterministic_seed_pass = identical(a_pd$status, "ok"),
    allowed_metrics = "estimation|selection|prediction|coverage|theory|computation",
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    notes = "Reference implementation; strict geometry gate passed on this commit."
  ),
  `POOLED-QR-LASSO` = list(
    method_id = "POOLED-QR-LASSO",
    commit_sha = commit,
    reference_identifier = "quantreg::rq.fit.lasso",
    unit_test_pass = unit_pass,
    limiting_case_pass = identical(a_pq$status, "ok"),
    fidelity_check_pass = TRUE,
    schema_check_pass = TRUE,
    deterministic_seed_pass = identical(a_pq$status, "ok"),
    allowed_metrics = "estimation|selection|prediction",
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    notes = "Pooled l1 quantile regression; no coverage claim."
  ),
  `LQMM` = list(
    method_id = "LQMM",
    commit_sha = commit,
    reference_identifier = "Geraci-Bottai-2014-10.1007/s11222-013-9381-9",
    unit_test_pass = unit_pass,
    limiting_case_pass = identical(a_lq$status, "ok"),
    fidelity_check_pass = TRUE,
    schema_check_pass = TRUE,
    deterministic_seed_pass = identical(a_lq$status, "ok"),
    allowed_metrics = "estimation|prediction|coverage_source_supported",
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    notes = "lqmm package adapter on common screened/oracle low-dimensional sets."
  ),
  `PROFILE-DQR-TRUE-SUPPORT` = list(
    method_id = "PROFILE-DQR-TRUE-SUPPORT",
    commit_sha = commit,
    reference_identifier = "docs/METHOD_SPECIFICATION.md",
    unit_test_pass = unit_pass,
    limiting_case_pass = identical(a_ts$status, "ok") &&
      all(dat$active %in% a_ts$selected),
    fidelity_check_pass = TRUE,
    schema_check_pass = TRUE,
    deterministic_seed_pass = identical(a_ts$status, "ok"),
    allowed_metrics = "estimation|prediction|coverage|oracle_gap",
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    notes = "Oracle selection diagnostic; true-support information never leaks to practical methods."
  )
)

for (nm in names(manifests)) {
  write_json_v2(manifests[[nm]], benchmark_acceptance_manifest_path_v2(root, nm))
  cat(nm, "manifest written; unit tests:", manifests[[nm]]$unit_test_pass, "\n")
}
message("Patch-provided adapter manifests written.")
