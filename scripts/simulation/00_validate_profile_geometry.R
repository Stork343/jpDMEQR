#!/usr/bin/env Rscript

full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "scripts/simulation/00_validate_profile_geometry.R"
root <- if (file.exists("R/profile_v2.R")) "." else normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", envir = environment())
source_v2_module(root, "simulation_v2", envir = environment())

out_dir <- file.path(root, "results", "logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rows <- list()
count <- 0L
for (seed in c(101L, 303L, 707L)) {
  for (tau in c(0.25, 0.5, 0.75)) {
    for (q in c(1L, 2L)) {
      for (h_mult in c(0.75, 1.25)) {
        dat <- generate_profile_qr_data_v2(
          n_clusters = 10, p = 4, s = 2, tau = tau, q = q,
          m_values = c(3, 4, 6), error_dist = if (seed == 707L) "t3" else "gaussian",
          seed = seed + 100L * q
        )
        beta_probe <- dat$beta0 + c(0.08, -0.04, 0.03, -0.02)
        Lambda <- if (q == 1L) 1.3 else matrix(c(1.4, 0.18, 0.18, 0.9), 2, 2)
        count <- count + 1L
        rows[[count]] <- validate_profile_geometry_v2(
          y = dat$y, X = dat$X, Z = dat$Z, cluster_id = dat$cluster_id,
          beta = beta_probe, tau = tau,
          h = h_mult * dat$n_clusters^(-1 / 3),
          lambda_gamma = Lambda, eps = 1e-5
        )
        rows[[count]]$seed <- seed
        rows[[count]]$h_multiplier <- h_mult
      }
    }
  }
}

res <- do.call(rbind, rows)
path <- file.path(out_dir, "profile_geometry_validation.csv")
write.csv(res, path, row.names = FALSE)
print(res)

thresholds <- list(
  gradient = 2e-4,
  hessian = 3e-3,
  identity = 1e-8,
  nuisance = 1e-7
)
failed <- with(res,
  gradient_max_error > thresholds$gradient |
  hessian_max_error > thresholds$hessian |
  profile_identity_error > thresholds$identity |
  max_nuisance_gradient > thresholds$nuisance |
  !nuisance_converged
)
if (any(failed)) {
  print(res[failed, , drop = FALSE])
  stop("Profile geometry validation failed. See ", path)
}

if (requireNamespace("CVXR", quietly = TRUE)) {
  dat <- generate_profile_qr_data_v2(15, 4, 2, 0.5, q = 1, seed = 909)
  comp <- profile_components_v2(
    dat$y, dat$X, dat$Z, dat$cluster_id, dat$beta0,
    tau = 0.5, h = 15^(-1 / 3), lambda_gamma = 1
  )
  d <- solve_dantzig_row_v2(comp$hessian, 1, mu_grid = c(0.05, 0.1, 0.2, 0.4))
  if (!d$feasible || d$residual > d$mu + 1e-6) stop("Dantzig geometry check failed.")
  write.csv(d$attempts, file.path(out_dir, "dantzig_geometry_validation.csv"), row.names = FALSE)
} else {
  warning("CVXR is not installed; Dantzig validation was skipped.")
}

message("All available profile-geometry checks passed.")
