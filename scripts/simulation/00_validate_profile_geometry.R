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

# 装配模式的对照：逐簇累加与整块装配的求和次序不同，整 Hessian、四份折 Hessian
# 与两式恒等矩阵都要落在冻结的恒等水平内
# (docs/METHOD_SPECIFICATION.md 第 4.4 节)。环境选中的模式由身份文件记录，
# 这里显式按两种模式各算一次。
assembly_cases <- list(
  list(n = 10L, p = 4L, s = 2L, tau = 0.5, q = 1L, m_values = 4L,
       m_rule = "uniform", error_dist = "gaussian", seed = 909L, lambda = 1.3),
  list(n = 10L, p = 4L, s = 2L, tau = 0.25, q = 2L, m_values = 2:8,
       m_rule = "geometric_imbalance", error_dist = "t3", seed = 1010L,
       lambda = matrix(c(1.4, 0.18, 0.18, 0.9), 2, 2)),
  list(n = 24L, p = 40L, s = 5L, tau = 0.5, q = 1L, m_values = 3:6,
       m_rule = "uniform", error_dist = "gaussian", seed = 1011L, lambda = 1.3),
  list(n = 24L, p = 40L, s = 5L, tau = 0.75, q = 2L, m_values = 2:8,
       m_rule = "geometric_imbalance", error_dist = "gaussian", seed = 1012L,
       lambda = matrix(c(2.2, -0.3, -0.3, 0.8), 2, 2)),
  list(n = 40L, p = 80L, s = 5L, tau = 0.5, q = 2L, m_values = 4L,
       m_rule = "uniform", error_dist = "gaussian", seed = 1013L, lambda = 1.3),
  list(n = 60L, p = 120L, s = 5L, tau = 0.25, q = 1L, m_values = 2:8,
       m_rule = "geometric_imbalance", error_dist = "t3", seed = 1014L,
       lambda = 1.3)
)
perturbation <- c(0.08, -0.04, 0.03, -0.02)
assembly_rows <- list()
for (ii in seq_along(assembly_cases)) {
  cs <- assembly_cases[[ii]]
  dat_a <- generate_profile_qr_data_v2(
    n_clusters = cs$n, p = cs$p, s = cs$s, tau = cs$tau, q = cs$q,
    m_values = cs$m_values, m_rule = cs$m_rule, error_dist = cs$error_dist,
    seed = cs$seed
  )
  beta_a <- dat_a$beta0
  beta_a[seq_along(perturbation)] <- beta_a[seq_along(perturbation)] + perturbation
  comp <- lapply(c("sequential", "batch"), function(mode) {
    profile_components_v2(
      dat_a$y, dat_a$X, dat_a$Z, dat_a$cluster_id, beta_a, cs$tau,
      cs$n^(-3 / 10), lambda_gamma = cs$lambda, need_hessian = TRUE,
      identity_check = TRUE, need_fold_sums = TRUE, assembly = mode
    )
  })
  hessian_difference <- max(abs(comp[[1]]$hessian - comp[[2]]$hessian))
  fold_difference <- max(vapply(seq_len(4L), function(f) {
    max(abs(comp[[1]]$hessian_fold_sums[[f]] - comp[[2]]$hessian_fold_sums[[f]]))
  }, numeric(1)))
  identity_difference <- max(abs(comp[[1]]$hessian_identity - comp[[2]]$hessian_identity))
  identity_error <- max(comp[[1]]$profile_identity_error, comp[[2]]$profile_identity_error)
  assembly_rows[[ii]] <- data.frame(
    case = ii, n_clusters = cs$n, p = cs$p, q = cs$q, tau = cs$tau,
    m_rule = cs$m_rule, error_dist = cs$error_dist,
    hessian_max_difference = hessian_difference,
    fold_hessian_max_difference = fold_difference,
    identity_matrix_max_difference = identity_difference,
    identity_error_max = identity_error,
    pass = max(hessian_difference, fold_difference, identity_difference,
               identity_error) <= thresholds$identity,
    stringsAsFactors = FALSE
  )
}
assembly <- do.call(rbind, assembly_rows)
write_atomic_csv_v2(assembly, file.path(out_dir, "assembly_mode_validation.csv"))
assembly_pass <- all(assembly$pass)

pass <- !any(failed) && dantzig_pass && assembly_pass
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
  assembly_modes_checked = TRUE,
  assembly_pass = assembly_pass,
  assembly_case_count = nrow(assembly),
  assembly_cases = assembly,
  thresholds = thresholds,
  solver_configuration = solver_configuration_v2(),
  maxima = list(
    gradient = max(res$gradient_max_error),
    hessian = max(res$hessian_max_error),
    identity = max(res$profile_identity_error),
    nuisance = max(res$max_nuisance_gradient),
    assembly_hessian = max(assembly$hessian_max_difference),
    assembly_fold_hessian = max(assembly$fold_hessian_max_difference),
    assembly_identity_matrix = max(assembly$identity_matrix_max_difference)
  ),
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
)
write_json_v2(manifest, geometry_manifest_path_v2(root))

if (any(failed)) {
  print(res[failed, , drop = FALSE])
  stop("Profile geometry validation failed. See ", out_dir)
}
if (!dantzig_pass) stop("Dantzig geometry validation failed.")
if (!assembly_pass) {
  print(assembly[!assembly$pass, , drop = FALSE])
  stop("Assembly-mode geometry validation failed. See ", out_dir)
}
message("Profile-geometry gate passed (strict=", strict, ").")
