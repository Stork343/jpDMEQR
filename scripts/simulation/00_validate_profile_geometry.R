#!/usr/bin/env Rscript

full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else
  "scripts/simulation/00_validate_profile_geometry.R"
root <- if (file.exists("R/profile_v2.R")) "." else
  normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", envir = environment())
source_v2_module(root, "simulation_v2", envir = environment())
source_v2_module(root, "metrics_v2", envir = environment())
source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"))

cli <- parse_cli_args_v2(commandArgs(trailingOnly = TRUE))
strict <- as_bool_cli_v2(cli$strict, FALSE)
seed_count <- as.integer(cli$`seed-count` %||% if (strict) 20L else 3L)
if (strict && seed_count < 20L) stop("Strict geometry requires at least 20 seeds.")
seeds <- 101L + 607L * seq_len(seed_count)
out_dir <- file.path(root, "results", "preflight", "geometry")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

thresholds <- list(
  gradient = 2e-4,
  hessian = 3e-3,
  identity = 1e-8,
  nuisance = 1e-7,
  dantzig = 1e-6
)
rows <- list(); count <- 0L
for (seed in seeds) {
  for (tau in c(0.25, 0.5, 0.75)) {
    for (q in c(1L, 2L)) {
      for (h_mult in c(0.75, 1, 1.25)) {
        balance_types <- if (strict) c("balanced", "unbalanced") else
          if ((seed + q) %% 2L) "balanced" else "unbalanced"
        for (balance_type in balance_types) {
          m_values <- if (balance_type == "balanced") 4L else 2:8
          m_rule <- if (balance_type == "balanced") "uniform" else "geometric_imbalance"
          dat <- generate_profile_qr_data_v2(
            n_clusters = 10,
            p = 4,
            s = 2,
            tau = tau,
            q = q,
            m_values = m_values,
            m_rule = m_rule,
            error_dist = if (seed %% 3L == 0L) "t3" else "gaussian",
            seed = seed + 100L * q + round(1000 * h_mult)
          )
          beta_probe <- dat$beta0 + c(0.08, -0.04, 0.03, -0.02)
          Lambda <- if (q == 1L) 1.3 else matrix(c(1.4, 0.18, 0.18, 0.9), 2, 2)
          count <- count + 1L
          res <- validate_profile_geometry_v2(
            y = dat$y,
            X = dat$X,
            Z = dat$Z,
            cluster_id = dat$cluster_id,
            beta = beta_probe,
            tau = tau,
            h = h_mult * dat$n_clusters^(-3 / 10),
            lambda_gamma = Lambda,
            eps = 1e-5
          )
          res$seed <- seed
          res$h_multiplier <- h_mult
          res$balance_type <- balance_type
          res$finite_difference_eps <- 1e-5
          rows[[count]] <- res
        }
      }
    }
  }
}

# Step-size sensitivity is required in strict mode but need not duplicate the
# entire grid. Use one balanced and one unbalanced case for every seed.
if (strict) {
  for (seed in seeds) {
    for (balance_type in c("balanced", "unbalanced")) {
      dat <- generate_profile_qr_data_v2(
        n_clusters = 10, p = 4, s = 2, tau = 0.5, q = 2,
        m_values = if (balance_type == "balanced") 4L else 2:8,
        m_rule = if (balance_type == "balanced") "uniform" else "geometric_imbalance",
        error_dist = "gaussian", seed = seed + 9191L
      )
      count <- count + 1L
      res <- validate_profile_geometry_v2(
        y = dat$y, X = dat$X, Z = dat$Z, cluster_id = dat$cluster_id,
        beta = dat$beta0 + c(0.08, -0.04, 0.03, -0.02),
        tau = 0.5, h = dat$n_clusters^(-3 / 10),
        lambda_gamma = matrix(c(1.4, 0.18, 0.18, 0.9), 2, 2),
        eps = 5e-6
      )
      res$seed <- seed
      res$h_multiplier <- 1
      res$balance_type <- balance_type
      res$finite_difference_eps <- 5e-6
      rows[[count]] <- res
    }
  }
}

res <- do.call(rbind, rows)
# Acceptance follows the frozen thresholds in docs/METHOD_SPECIFICATION.md
# (Section 11): nuisance gradient maximum error below 1e-7. The solver-internal
# convergence flag is stricter (1e-8) and is kept as a diagnostic column rather
# than a gate criterion; small-cluster t3 cases can sit on a numerical plateau
# between 1e-8 and 1e-7 that still satisfies the frozen acceptance threshold.
failed <- with(res,
  gradient_max_error > thresholds$gradient |
    hessian_max_error > thresholds$hessian |
    profile_identity_error > thresholds$identity |
    max_nuisance_gradient > thresholds$nuisance
)
write_atomic_csv_v2(res, file.path(out_dir, "profile_geometry_validation.csv"))

if (!requireNamespace("CVXR", quietly = TRUE)) {
  stop("CVXR is mandatory for the geometry gate; install dependencies first.")
}
dat_d <- generate_profile_qr_data_v2(18, 5, 2, 0.5, q = 2, seed = 909L)
comp_d <- profile_components_v2(
  dat_d$y, dat_d$X, dat_d$Z, dat_d$cluster_id, dat_d$beta0,
  tau = 0.5, h = 18^(-3 / 10),
  lambda_gamma = matrix(c(1.4, 0.18, 0.18, 0.9), 2, 2)
)
dantzig <- solve_dantzig_row_v2(
  comp_d$hessian, coordinate = 1,
  mu_grid = c(0.05, 0.1, 0.2, 0.4, 0.8)
)
dantzig_pass <- isTRUE(dantzig$feasible) &&
  dantzig$residual <= dantzig$mu + thresholds$dantzig
write_atomic_csv_v2(dantzig$attempts, file.path(out_dir, "dantzig_geometry_validation.csv"))

pass <- !any(failed) && dantzig_pass
manifest <- list(
  pass = pass,
  strict = strict,
  commit_sha = current_commit_v2(root),
  branch = current_branch_v2(root),
  seed_count = length(seeds),
  case_count = nrow(res),
  tau_values = c(0.25, 0.5, 0.75),
  q_values = c(1, 2),
  h_multipliers = c(0.75, 1, 1.25),
  balance_types = unique(res$balance_type),
  finite_difference_steps = sort(unique(res$finite_difference_eps)),
  non_diagonal_lambda_checked = TRUE,
  dantzig_checked = TRUE,
  dantzig_pass = dantzig_pass,
  dantzig_residual = dantzig$residual,
  dantzig_mu = dantzig$mu,
  thresholds = thresholds,
  maxima = list(
    gradient = max(res$gradient_max_error),
    hessian = max(res$hessian_max_error),
    identity = max(res$profile_identity_error),
    nuisance = max(res$max_nuisance_gradient)
  ),
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
)
write_json_v2(manifest, geometry_manifest_path_v2(root))

if (any(failed)) {
  print(res[failed, , drop = FALSE])
  stop("Profile geometry validation failed. See ", out_dir)
}
if (!dantzig_pass) stop("Dantzig geometry validation failed.")
message("Profile-geometry gate passed (strict=", strict, ").")
