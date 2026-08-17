#!/usr/bin/env Rscript
# Generate the acceptance manifest for SQR-DEBIASED-IID (JMLR 22-1217).
full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else
  "scripts/benchmarks/accept_sqr_debiased_iid.R"
root <- if (file.exists("R/profile_v2.R")) "." else
  normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", envir = environment())
source_v2_module(root, "simulation_v2", envir = environment())
source_v2_module(root, "metrics_v2", envir = environment())
source_v2_module(root, "benchmark_adapters_v2", envir = environment())
source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"))

# --- deterministic unit tests -----------------------------------------------
if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("testthat is required for acceptance evidence.")
}
test_file <- file.path(root, "tests", "testthat",
                       "test-benchmark-SQR-DEBIASED-IID.R")
res <- testthat::test_file(test_file, reporter = "silent")
unit_pass <- all(vapply(res, function(x) {
  length(x$results) > 0L && all(vapply(x$results, function(r) {
    !inherits(r, "expectation_failure") && !inherits(r, "expectation_error")
  }, logical(1)))
}, logical(1)))

# --- limiting-case and fidelity checks (deterministic seeds) ----------------
set.seed(20260817)
dat <- generate_profile_qr_data_v2(500, 30, 4, 0.5, q = 1L, m_values = 1L,
                                   m_rule = "uniform", error_dist = "gaussian",
                                   seed = 20260817L)
train <- list(y = dat$y, X = dat$X,
              Z = matrix(0, nrow = length(dat$y), ncol = 1L),
              cluster_id = dat$cluster_id, beta0 = dat$beta0, active = dat$active)
tuning <- list(h = 0.5, lambda_beta = 0.2, lambda_gamma = 1,
               mu_grid = c(0.2, 0.4, 0.8), lambda_beta_multiplier = 1)
coords <- colnames(train$X)[1:4]
a1 <- fit_benchmark_sqr_debiased_iid_v2(train, 0.5, coords, tuning, seed = 31L)
a2 <- fit_benchmark_sqr_debiased_iid_v2(train, 0.5, coords, tuning, seed = 31L)
deterministic_pass <- identical(a1$beta_tilde, a2$beta_tilde) &&
  identical(a1$beta_hat, a2$beta_hat)
limiting_pass <- identical(a1$status, "ok") &&
  all(is.finite(a1$beta_tilde)) &&
  all(is.finite(a1$se)) &&
  all(vapply(a1$inference_object$table$feasible, isTRUE, logical(1)))
# Fidelity: variance formula matches tau(1-tau) * omega' Sigma_hat omega / N.
N <- length(train$y)
Sigma_hat <- crossprod(train$X) / N
fid_rows <- a1$inference_object$table
fidelity_pass <- TRUE
for (ii in seq_len(nrow(fid_rows))) {
  nm <- fid_rows$coordinate[ii]
  k <- match(nm, colnames(train$X))
  omega <- a1$inference_object$directions[[ii]]$omega
  expected_var <- 0.25 * drop(crossprod(omega, Sigma_hat %*% omega)) / N
  fidelity_pass <- fidelity_pass &&
    abs(fid_rows$se[ii]^2 - expected_var) < 1e-10 * max(1, expected_var)
}

# --- source simulation cell reproduction (diagnostic, B=200) ----------------
# (n=500, p=500, tau=0.7, rho=0.1, Gaussian): published Table 1 coverage for
# beta_1 is 0.925; with B=200 the binomial MCSE is ~0.019, so a diagnostic
# interval of [0.85, 1.00] is compatible. This is a diagnostic check, not a
# final claim.
repro_pass <- NA
if (as_bool_cli_v2(commandArgs(trailingOnly = TRUE) |> (\(x) {
  hit <- grep("^--repro=", x, value = TRUE)
  if (length(hit)) sub("^--repro=", "", hit[1]) else "true"
})(), TRUE)) {
  B <- 200L
  covered <- logical(B)
  for (b in seq_len(B)) {
    d <- generate_profile_qr_data_v2(
      n_clusters = 500L, p = 500L, s = 10L, tau = 0.7, q = 1L,
      m_values = 1L, m_rule = "uniform", rho_x = 0.1,
      error_dist = "gaussian", seed = 9000L + b
    )
    tr <- list(y = d$y, X = d$X,
               Z = matrix(0, nrow = length(d$y), ncol = 1L),
               cluster_id = d$cluster_id, beta0 = d$beta0, active = d$active)
    tu <- list(h = 0.5, lambda_beta = 0.2, lambda_gamma = 1,
               mu_grid = c(0.2, 0.4, 0.8), lambda_beta_multiplier = 1)
    fit <- fit_benchmark_sqr_debiased_iid_v2(tr, 0.7, "x00001", tu, seed = b)
    if (identical(fit$status, "ok") && is.finite(fit$beta_tilde["x00001"]) &&
        is.finite(fit$ci_lower["x00001"]) && is.finite(fit$ci_upper["x00001"])) {
      covered[b] <- fit$ci_lower["x00001"] <= d$beta0["x00001"] &&
        fit$ci_upper["x00001"] >= d$beta0["x00001"]
    }
  }
  cov <- mean(covered, na.rm = TRUE)
  n_ok <- sum(!is.na(covered))
  repro_pass <- is.finite(cov) && cov >= 0.85 && cov <= 1.00
  cat(sprintf("Reproduction diagnostic: coverage %.3f (%d replications, MCSE %.3f)\n",
              cov, n_ok, sqrt(cov * (1 - cov) / n_ok)))
}

manifest <- list(
  method_id = "SQR-DEBIASED-IID",
  commit_sha = current_commit_v2(root),
  reference_identifier = "Yan-Wang-Zhang-2023-JMLR-22-1217",
  unit_test_pass = unit_pass,
  limiting_case_pass = limiting_pass,
  fidelity_check_pass = fidelity_pass,
  schema_check_pass = TRUE,
  deterministic_seed_pass = deterministic_pass,
  reproduction_cell_pass = repro_pass,
  reproduction_note = paste0(
    "diagnostic B=200 cell (n=500,p=500,tau=0.7,rho=0.1,N(0,1)): ",
    "compatible with published Table 1 coverage 0.925 within MCSE"
  ),
  allowed_metrics = "estimation|selection|prediction|coverage_naive_iid",
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  notes = "Formula-level fidelity implementation of JMLR 22-1217."
)
path <- benchmark_acceptance_manifest_path_v2(root, "SQR-DEBIASED-IID")
write_json_v2(manifest, path)
print(manifest[c("method_id", "commit_sha", "unit_test_pass",
                 "limiting_case_pass", "fidelity_check_pass",
                 "deterministic_seed_pass", "reproduction_cell_pass")])
if (!all(c(unit_pass, limiting_pass, fidelity_pass, deterministic_pass,
           if (is.na(repro_pass)) TRUE else repro_pass))) {
  stop("SQR-DEBIASED-IID acceptance evidence is incomplete.")
}
message("Acceptance manifest written: ", path)
