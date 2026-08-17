# Fidelity implementation of BIAS-ADJ-LQMM.
#
# Source method: Battagliola, M., Sørensen, H., Tolver, A., Staicu, A.-M.
# (2022). "A bias adjusted estimator in quantile regression for clustered
# data". Econometrics and Statistics, DOI 10.1016/j.ecosta.2021.07.003.
# arXiv:2202.11501. Official R code: bias_adj.zip (helle.sites.ku.dk).
#
# The method is the two-step estimator with RW (resampling + wild) bootstrap
# bias correction:
#   Step 1: LQMM (Geraci & Bottai) random-intercept fit; take the centred BLP
#           random intercepts tilde{u}_i.
#   Step 2: standard quantile regression of y - Z' tilde{u}_i on X gives the
#           two-step estimate beta_hat.
#   Bootstrap (B replications):
#     w_ij ~ {2(1-tau) w.p. 1-tau; -2tau w.p. tau} (wild, zero tau-quantile);
#     u_i* ~ resample from {tilde{u}_i} (cluster-level, with replacement);
#     Y*_ij = X_ij' beta_hat + Z_ij' u_i* + w_ij |residual_ij|;
#     oracle replicate: rq of Y* - Z' u_i* on X;
#     two-step replicate: re-run Step 1 (LQMM) + Step 2 (rq) on Y*.
#   Adjusted estimate: beta_adj = 2 beta_hat - mean(beta_two_step*).
#   Adjusted SE: SE_adj = SE_obs * SD(two-step*) / SD(oracle*).
#   SE-adjusted CI: beta_adj +/- z_{1-alpha/2} SE_adj.
#   Basic bootstrap CI: (2 beta_hat - q_{1-alpha/2}, 2 beta_hat - q_{alpha/2}).

fit_benchmark_bias_adj_lqmm_v2 <- function(train,
                                           tau,
                                           target_coords,
                                           tuning,
                                           seed,
                                           context = list(),
                                           control = list()) {
  if (!requireNamespace("lqmm", quietly = TRUE) ||
      !requireNamespace("quantreg", quietly = TRUE)) {
    return(benchmark_failure_v2("BIAS-ADJ-LQMM", "dependency",
                                "lqmm and quantreg are required."))
  }
  set.seed(seed)
  start <- proc.time()[[3]]
  X <- as.matrix(train$X)
  y <- as.numeric(train$y)
  Z <- as.matrix(train$Z)
  cluster_id <- as.character(train$cluster_id)
  if (nrow(X) != length(y) || nrow(Z) != length(y) ||
      length(cluster_id) != length(y)) {
    return(benchmark_failure_v2("BIAS-ADJ-LQMM", "schema",
                                "Dimension mismatch in training data."))
  }
  if (ncol(Z) > 1L) {
    # The source method (ran_int_bias_adj.R) fits a random intercept only
    # (random = ~1), regardless of the generated nuisance dimension. Keep
    # the first column of Z (the intercept direction) as the random-effect
    # design; this is the faithful source behaviour rather than a failure.
    Z <- Z[, 1L, drop = FALSE]
  }
  B <- as.integer(control$B %||% 100L)
  nK <- as.integer(control$nK %||% 15L)
  level <- control$ci_level %||% 0.95
  alpha <- 1 - level
  zcrit <- stats::qnorm(1 - alpha / 2)

  df <- data.frame(
    y = y,
    cluster = factor(cluster_id),
    X,
    check.names = FALSE
  )
  fixed_names <- colnames(X)
  # Random-intercept LQMM formula: y ~ X (no fixed intercept is forced unless
  # the design contains one; the source method always includes the intercept
  # via lqmm's default formula handling).
  fixed_formula <- stats::as.formula(
    paste("y ~ 0 +", paste(fixed_names, collapse = " + "))
  )

  # Step 1: LQMM with Gauss-Hermite quadrature (source defaults nK=15,
  # type="normal", control method="df").
  step1 <- tryCatch(
    lqmm::lqmm(
      fixed = fixed_formula,
      random = ~1,
      group = cluster,
      tau = tau,
      data = df,
      nK = nK,
      type = control$lqmm_type %||% "normal",
      control = lqmm::lqmmControl(
        method = "df", LP_max_iter = 1000, verbose = FALSE
      )
    ),
    error = function(e) e
  )
  if (inherits(step1, "error")) {
    return(benchmark_failure_v2("BIAS-ADJ-LQMM", "nuisance_fit",
                                conditionMessage(step1),
                                proc.time()[[3]] - start))
  }

  # Centred BLP random intercepts per cluster.
  ran <- tryCatch(lqmm::ranef(step1), error = function(e) NULL)
  if (is.null(ran) || !"(Intercept)" %in% colnames(ran)) {
    return(benchmark_failure_v2("BIAS-ADJ-LQMM", "nuisance_fit",
                                "Could not extract lqmm random effects.",
                                proc.time()[[3]] - start))
  }
  u_est_cluster <- as.numeric(ran[, "(Intercept)"]) -
    mean(as.numeric(ran[, "(Intercept)"]))
  cluster_levels <- levels(df$cluster)
  if (length(u_est_cluster) != length(cluster_levels)) {
    return(benchmark_failure_v2("BIAS-ADJ-LQMM", "nuisance_fit",
                                "Random-effect/cluster alignment failed.",
                                proc.time()[[3]] - start))
  }
  u_est <- u_est_cluster[match(df$cluster, cluster_levels)]

  # Step 2: standard QR on centred response.
  df_step2 <- df
  df_step2$y <- y - as.numeric(Z[, 1] * u_est)
  step2 <- tryCatch(
    quantreg::rq(y ~ 0 + ., data = df_step2, tau = tau),
    error = function(e) e
  )
  if (inherits(step2, "error")) {
    return(benchmark_failure_v2("BIAS-ADJ-LQMM", "penalised_fit",
                                conditionMessage(step2),
                                proc.time()[[3]] - start))
  }
  beta_two_step <- stats::coef(step2)
  if (is.null(names(beta_two_step)) ||
      !all(fixed_names %in% names(beta_two_step))) {
    return(benchmark_failure_v2("BIAS-ADJ-LQMM", "coefficient_extraction",
                                "rq coefficient names are misaligned.",
                                proc.time()[[3]] - start))
  }
  beta_two_step <- beta_two_step[fixed_names]
  se_obs <- tryCatch(
    stats::coef(summary(step2, se = control$se_method %||% "nid"))[, 2],
    error = function(e) rep(NA_real_, length(fixed_names))
  )
  se_obs <- se_obs[fixed_names]

  # Residuals from the two-step fit (with centred response).
  res <- as.numeric(stats::residuals(step2))
  fitted_val <- as.numeric(stats::fitted(step2))

  # Bootstrap.
  boot_oracle <- matrix(NA_real_, B, length(fixed_names))
  boot_two_step <- matrix(NA_real_, B, length(fixed_names))
  colnames(boot_oracle) <- colnames(boot_two_step) <- fixed_names
  N <- length(cluster_levels)
  nobs <- length(y)
  ok_boot <- 0L

  for (b in seq_len(B)) {
    w <- sample(
      c(2 * (1 - tau), -2 * tau), size = nobs, replace = TRUE,
      prob = c(1 - tau, tau)
    )
    u_boot_cluster <- sample(u_est_cluster, size = N, replace = TRUE)
    u_boot <- u_boot_cluster[match(df$cluster, cluster_levels)]
    y_boot <- fitted_val + as.numeric(Z[, 1] * u_boot) + w * abs(res)

    # Oracle replicate: rq of y_boot - Z' u_boot on X.
    df_oracle <- df
    df_oracle$y <- y_boot - as.numeric(Z[, 1] * u_boot)
    fit_oracle <- tryCatch(
      quantreg::rq(y ~ 0 + ., data = df_oracle, tau = tau),
      error = function(e) NULL
    )
    if (is.null(fit_oracle)) next
    co <- stats::coef(fit_oracle)
    if (is.null(names(co)) || !all(fixed_names %in% names(co))) next
    boot_oracle[b, ] <- co[fixed_names]

    # Two-step replicate: re-run LQMM then rq on y_boot.
    df_boot <- df
    df_boot$y <- y_boot
    step1b <- tryCatch(
      lqmm::lqmm(
        fixed = fixed_formula, random = ~1, group = cluster, tau = tau,
        data = df_boot, nK = nK,
        type = control$lqmm_type %||% "normal",
        control = lqmm::lqmmControl(
          method = "df", LP_max_iter = 1000, verbose = FALSE
        )
      ),
      error = function(e) NULL
    )
    if (is.null(step1b)) next
    ranb <- tryCatch(lqmm::ranef(step1b), error = function(e) NULL)
    if (is.null(ranb) || !"(Intercept)" %in% colnames(ranb)) next
    u_estb_cluster <- as.numeric(ranb[, "(Intercept)"]) -
      mean(as.numeric(ranb[, "(Intercept)"]))
    if (length(u_estb_cluster) != N) next
    u_estb <- u_estb_cluster[match(df_boot$cluster, cluster_levels)]
    df_step2b <- df_boot
    df_step2b$y <- y_boot - as.numeric(Z[, 1] * u_estb)
    step2b <- tryCatch(
      quantreg::rq(y ~ 0 + ., data = df_step2b, tau = tau),
      error = function(e) NULL
    )
    if (is.null(step2b)) next
    cob <- stats::coef(step2b)
    if (is.null(names(cob)) || !all(fixed_names %in% names(cob))) next
    boot_two_step[b, ] <- cob[fixed_names]
    ok_boot <- ok_boot + 1L
  }

  min_ok <- control$min_bootstrap_fraction %||% 0.5
  if (ok_boot < B * min_ok) {
    return(benchmark_failure_v2(
      "BIAS-ADJ-LQMM", "resampling",
      sprintf("Only %d/%d bootstrap replications succeeded.", ok_boot, B),
      proc.time()[[3]] - start
    ))
  }

  # Bias adjustment (Efron & Tibshirani bootstrap bias correction).
  beta_adj <- 2 * beta_two_step - colMeans(boot_two_step, na.rm = TRUE)

  # SE adjustment: SE_adj = SE_obs * SD(two-step*)/SD(oracle*).
  sd_two <- apply(boot_two_step, 2, stats::sd, na.rm = TRUE)
  sd_oracle <- apply(boot_oracle, 2, stats::sd, na.rm = TRUE)
  ratio <- sd_two / sd_oracle
  ratio[!is.finite(ratio) | ratio <= 0] <- 1
  se_adj <- se_obs * ratio
  se_adj[!is.finite(se_adj)] <- NA_real_

  # SE-adjusted CI (recommended) and basic bootstrap CI.
  ci_lower <- beta_adj - zcrit * se_adj
  ci_upper <- beta_adj + zcrit * se_adj
  q_lo <- apply(boot_two_step, 2, stats::quantile, probs = alpha / 2, na.rm = TRUE)
  q_hi <- apply(boot_two_step, 2, stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE)
  basic_lower <- 2 * beta_two_step - q_hi
  basic_upper <- 2 * beta_two_step - q_lo

  elapsed <- proc.time()[[3]] - start
  # Report inference only for the requested target coordinates (the full
  # coefficient vector is still returned in beta_hat).
  requested <- if (is.character(target_coords)) {
    intersect(as.character(target_coords), fixed_names)
  } else {
    idx <- as.integer(target_coords)
    idx <- idx[idx >= 1L & idx <= length(fixed_names)]
    fixed_names[idx]
  }
  if (!length(requested)) {
    return(benchmark_failure_v2("BIAS-ADJ-LQMM", "schema",
                                "No valid target coordinates.",
                                proc.time()[[3]] - start))
  }
  keep <- fixed_names %in% requested
  inf_tab <- data.frame(
    coordinate = fixed_names[keep],
    index = which(keep),
    feasible = is.finite(se_adj[keep]),
    beta_hat = as.numeric(beta_two_step[keep]),
    beta_tilde = as.numeric(beta_adj[keep]),
    se = as.numeric(se_adj[keep]),
    ci_lower = as.numeric(ci_lower[keep]),
    ci_upper = as.numeric(ci_upper[keep]),
    mu = NA_real_,
    dantzig_residual = NA_real_,
    omega_l1 = NA_real_,
    omega_l2 = NA_real_,
    basic_ci_lower = as.numeric(basic_lower[keep]),
    basic_ci_upper = as.numeric(basic_upper[keep]),
    stringsAsFactors = FALSE
  )
  ans <- list(
    method_id = "BIAS-ADJ-LQMM",
    status = if (all(is.finite(beta_adj))) "ok" else "warning",
    failure_stage = "",
    failure_message = "",
    beta_hat = beta_two_step,
    beta_tilde = beta_adj[keep],
    se = se_adj[keep],
    ci_lower = ci_lower[keep],
    ci_upper = ci_upper[keep],
    selected = which(abs(beta_two_step) > (control$selection_tol %||% 1e-8)),
    fit_object = list(
      two_step = step2,
      lqmm = step1,
      u_blp = u_est_cluster,
      bootstrap_success = ok_boot,
      B = B,
      beta_adj = beta_adj,
      se_obs = se_obs,
      sd_two_step = sd_two,
      sd_oracle = sd_oracle
    ),
    inference_object = list(
      table = inf_tab,
      ci_level = level,
      interval_type = "se_adjusted_and_basic_bootstrap"
    ),
    runtime_sec = elapsed,
    converged = ok_boot >= B * 0.9,
    kkt_residual = NA_real_,
    warning_messages = if (ok_boot < B) {
      sprintf("Bootstrap completion %d/%d.", ok_boot, B)
    } else character(),
    implementation_version = "battagliola-rw-bootstrap",
    target_scope = "source_defined_clustered_quantile"
  )
  benchmark_add_metadata_v2(
    ans,
    reference_identifier = "Battagliola-et-al-2022-10.1016/j.ecosta.2021.07.003",
    fidelity_status = "formula_level_implemented"
  )
}
