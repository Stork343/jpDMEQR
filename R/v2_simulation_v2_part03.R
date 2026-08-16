approximate_profile_target_v2 <- function(config,
                                          target_columns,
                                          n_population = 100000L,
                                          seed = 87001L,
                                          fit_control = list(max_iter = 500L,
                                                             beta_tol = 1e-8,
                                                             kkt_tol = 1e-7)) {
  if (!exists("fit_profile_lasso_v2", mode = "function")) {
    stop("Source R/profile_v2.R before approximating the profile target.")
  }
  cfg <- config
  cfg$n_clusters <- as.integer(n_population)
  cfg$p <- max(target_columns)
  dat <- simulate_from_config_v2(cfg, seed = seed)
  if (identical(cfg$fit_random_effects %||% "", "intercept") && ncol(dat$Z) > 1L) {
    dat$Z <- dat$Z[, 1L, drop = FALSE]
  }
  X_target <- dat$X[, target_columns, drop = FALSE]
  fit <- fit_profile_lasso_v2(
    y = dat$y, X = X_target, Z = dat$Z, cluster_id = dat$cluster_id,
    tau = cfg$tau, h = n_population^(-1 / 3),
    lambda_beta = 0, lambda_gamma = cfg$lambda_gamma,
    penalty_factor = rep(0, length(target_columns)),
    control = fit_control
  )
  list(
    beta_star_mc = fit$beta,
    converged = fit$converged,
    kkt_residual = fit$kkt_residual,
    n_population = n_population,
    seed = seed,
    target_columns = target_columns,
    fit = fit
  )
}
