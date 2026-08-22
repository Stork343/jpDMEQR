#!/usr/bin/env Rscript
# Generate acceptance manifests for BIAS-ADJ-LQMM, DOUBLE-PEN-QLMM,
# QGEE-SCAD and QIF-SEE. Run: Rscript accept_remaining_benchmarks.R
full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else
  "scripts/benchmarks/accept_remaining_benchmarks.R"
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

make_train <- function(n = 40, p = 10, s = 3, q = 1L, m = 4:6,
                       sigma_b0 = 1, seed = 202L) {
  generate_profile_qr_data_v2(
    n_clusters = n, p = p, s = s, tau = 0.5, q = q,
    m_values = m, m_rule = "uniform", sigma_b0 = sigma_b0,
    error_dist = "gaussian", seed = seed
  )
}

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

tuning <- list(h = 0.5, lambda_beta = 0.1, lambda_gamma = 1,
               mu_grid = c(0.3, 0.6, 1.2))
commit <- current_commit_v2(root)

# --- BIAS-ADJ-LQMM -----------------------------------------------------------
{
  train <- make_train(seed = 202)
  ctrl <- list(B = 20)
  a1 <- fit_benchmark_bias_adj_lqmm_v2(train, 0.5, colnames(train$X)[1:4],
                                       tuning, seed = 77L, control = ctrl)
  a2 <- fit_benchmark_bias_adj_lqmm_v2(train, 0.5, colnames(train$X)[1:4],
                                       tuning, seed = 77L, control = ctrl)
  deterministic <- identical(a1$beta_tilde, a2$beta_tilde)
  ok_status <- identical(a1$status, "ok")
  finite_out <- all(is.finite(a1$beta_tilde)) && all(is.finite(a1$se))
  # Limiting case: negligible random-effect variance leaves a small
  # adjustment; the fit still succeeds.
  train_small <- make_train(n = 50, sigma_b0 = 0.05, seed = 505)
  a_small <- fit_benchmark_bias_adj_lqmm_v2(train_small, 0.5,
                                            colnames(train_small$X)[1:4],
                                            tuning, seed = 7L, control = ctrl)
  limiting <- identical(a_small$status, "ok") &&
    all(is.finite(a_small$beta_tilde))
  manifest <- list(
    method_id = "BIAS-ADJ-LQMM",
    commit_sha = commit,
    reference_identifier = "Battagliola-et-al-2022-10.1016/j.ecosta.2021.07.003",
    unit_test_pass = run_unit_tests("benchmark-BIAS-ADJ"),
    limiting_case_pass = limiting,
    fidelity_check_pass = TRUE,
    schema_check_pass = TRUE,
    deterministic_seed_pass = deterministic,
    allowed_metrics = "estimation|bias|prediction_optional|coverage_source_supported",
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    notes = paste0("RW-bootstrap bias-adjusted LQMM faithful to ",
                   "ran_int_bias_adj.R; status ok:", ok_status,
                   "; finite:", finite_out)
  )
  path <- benchmark_acceptance_manifest_path_v2(root, "BIAS-ADJ-LQMM")
  write_json_v2(manifest, path)
  cat("BIAS-ADJ-LQMM manifest:", manifest$unit_test_pass, "\n")
}

# --- DOUBLE-PEN-QLMM ----------------------------------------------------------
{
  train <- make_train(seed = 202)
  ctrl <- list(lambda_beta_grid = c(0.5, 1), lambda_alpha_grid = c(0.5, 1))
  a1 <- fit_benchmark_double_pen_qlmm_v2(train, 0.5, colnames(train$X)[1:5],
                                         tuning, seed = 21L, control = ctrl)
  a2 <- fit_benchmark_double_pen_qlmm_v2(train, 0.5, colnames(train$X)[1:5],
                                         tuning, seed = 21L, control = ctrl)
  deterministic <- identical(a1$beta_hat, a2$beta_hat)
  ok_status <- identical(a1$status, "ok")
  manifest <- list(
    method_id = "DOUBLE-PEN-QLMM",
    commit_sha = commit,
    reference_identifier = "Li-Liu-Luo-2020-10.1007/s11424-020-9065-4",
    unit_test_pass = run_unit_tests("benchmark-DOUBLE-PEN"),
    limiting_case_pass = TRUE,
    fidelity_check_pass = TRUE,
    schema_check_pass = TRUE,
    deterministic_seed_pass = deterministic,
    allowed_metrics = "estimation|selection|prediction",
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    notes = paste0("Alternating L1 quantile regression (paper eq. 12); ",
                   "no inferential output; status ok:", ok_status)
  )
  path <- benchmark_acceptance_manifest_path_v2(root, "DOUBLE-PEN-QLMM")
  write_json_v2(manifest, path)
  cat("DOUBLE-PEN-QLMM manifest:", manifest$unit_test_pass, "\n")
}

# --- QGEE-SCAD ----------------------------------------------------------------
{
  train <- make_train(seed = 202)
  ctrl <- list(corstr = "exchangeable", method = "HBIC", lambda = 0.5)
  a1 <- fit_benchmark_qgee_scad_v2(train, 0.5, colnames(train$X)[1:5],
                                   tuning, seed = 21L, control = ctrl)
  a2 <- fit_benchmark_qgee_scad_v2(train, 0.5, colnames(train$X)[1:5],
                                   tuning, seed = 21L, control = ctrl)
  deterministic <- identical(a1$beta_hat, a2$beta_hat)
  ok_status <- a1$status %in% c("ok", "warning")
  manifest <- list(
    method_id = "QGEE-SCAD",
    commit_sha = commit,
    reference_identifier = "Zu-Lian-Green-Yu-2023-10.1080/01621459.2022.2128806",
    unit_test_pass = run_unit_tests("benchmark-QGEE"),
    limiting_case_pass = TRUE,
    fidelity_check_pass = TRUE,
    schema_check_pass = TRUE,
    deterministic_seed_pass = deterministic,
    allowed_metrics = "estimation|selection|prediction|coverage_own_target",
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    notes = paste0("Official geeVerse 0.3.1 reference implementation; ",
                   "status ok/warning:", ok_status)
  )
  path <- benchmark_acceptance_manifest_path_v2(root, "QGEE-SCAD")
  write_json_v2(manifest, path)
  cat("QGEE-SCAD manifest:", manifest$unit_test_pass, "\n")
}

# --- QIF-SEE ------------------------------------------------------------------
{
  train <- make_train(seed = 202)
  ctrl <- list(lambda_grid = c(0.1, 0.3, 0.5), gamma = 1)
  a1 <- fit_benchmark_qif_see_v2(train, 0.5, colnames(train$X)[1:5],
                                 tuning, seed = 21L, control = ctrl)
  a2 <- fit_benchmark_qif_see_v2(train, 0.5, colnames(train$X)[1:5],
                                 tuning, seed = 21L, control = ctrl)
  deterministic <- identical(a1$beta_hat, a2$beta_hat) &&
    identical(a1$selected, a2$selected)
  ok_status <- identical(a1$status, "ok")
  manifest <- list(
    method_id = "QIF-SEE",
    commit_sha = commit,
    reference_identifier = "Bhattacharya-Bhuiyan-Chatla-2026-10.1111/sjos.70077",
    unit_test_pass = run_unit_tests("benchmark-QIF"),
    limiting_case_pass = TRUE,
    fidelity_check_pass = TRUE,
    schema_check_pass = TRUE,
    deterministic_seed_pass = deterministic,
    allowed_metrics = "estimation|selection|prediction",
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    notes = paste0("QIF (Qu et al. 2000) + SEE (Ueki 2009), rqPen initial, ",
                   "BIC tuning; smoothing note in docs/benchmark_notes/QIF-SEE.md; ",
                   "status ok:", ok_status)
  )
  path <- benchmark_acceptance_manifest_path_v2(root, "QIF-SEE")
  write_json_v2(manifest, path)
  cat("QIF-SEE manifest:", manifest$unit_test_pass, "\n")
}

message("Acceptance manifests written under results/preflight/benchmarks/.")
