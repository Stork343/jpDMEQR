fit_benchmark_lqmm_v2 <- function(train, tau, target_coords, tuning, seed,
                                  random_slope = FALSE, control = list()) {
  if (!requireNamespace("lqmm", quietly = TRUE)) {
    return(benchmark_failure_v2("LQMM", "dependency", "lqmm is not installed"))
  }
  set.seed(seed)
  X <- as.matrix(train$X)
  if (ncol(X) > (control$max_features %||% 30L)) {
    return(benchmark_failure_v2(
      "LQMM", "dimension",
      "LQMM adapter requires a common screened or oracle set of at most max_features."
    ))
  }
  df <- data.frame(
    y = train$y,
    cluster = factor(train$cluster_id),
    time = if (ncol(train$Z) >= 2) train$Z[, 2] else 0,
    X,
    check.names = FALSE
  )
  # Standalone grouping vector: lqmm re-evaluates `group` in the formula
  # environment, so passing df[["cluster"]] collides with the data.frame
  # symbol `df`. A distinct variable name avoids the closure conflict.

  fixed_formula <- stats::as.formula(
    paste("y ~ 0 +", paste(colnames(X), collapse = " + "))
  )
  random_formula <- if (random_slope) ~ 1 + time else ~ 1
  start <- proc.time()[[3]]
  fit <- tryCatch(
    lqmm::lqmm(
      fixed = fixed_formula,
      random = random_formula,
      group = "cluster",
      tau = tau,
      data = df,
      covariance = "pdSymm",
      nK = control$nK %||% 7L,
      control = lqmm::lqmmControl(LP_tol_ll = control$tol %||% 1e-6)
    ),
    error = function(e) e
  )
  elapsed <- proc.time()[[3]] - start
  if (inherits(fit, "error")) {
    return(benchmark_failure_v2("LQMM", "penalised_fit", conditionMessage(fit), elapsed))
  }
  # lqmm >= 1.5.8 removed the fixef export; coef() returns the fixed effects.
  beta <- tryCatch(as.numeric(stats::coef(fit)$fixed), error = function(e) NULL)
  if (is.null(beta)) beta <- tryCatch(as.numeric(stats::coef(fit)), error = function(e) NULL)
  if (is.null(beta)) {
    return(benchmark_failure_v2("LQMM", "coefficient_extraction",
                                "Could not extract fixed effects.", elapsed))
  }
  names(beta) <- colnames(X)[seq_along(beta)]
  list(
    method_id = "LQMM",
    status = "ok",
    failure_stage = "",
    failure_message = "",
    beta_hat = beta,
    beta_tilde = NULL,
    se = NULL,
    selected = which(abs(beta) > (control$selection_tol %||% 1e-8)),
    fit_object = fit,
    runtime_sec = elapsed,
    converged = TRUE,
    kkt_residual = NA_real_,
    warning_messages = character(),
    implementation_version = "lqmm"
  )
}

fit_benchmark_by_id_v2 <- function(method_id, train, tau, target_coords,
                                   tuning, seed, context = list(), control = list()) {
  switch(
    method_id,
    `PROFILE-DQR` = fit_benchmark_profile_dqr_v2(train, tau, target_coords, tuning, seed, control),
    `POOLED-QR-LASSO` = fit_benchmark_pooled_qr_lasso_v2(train, tau, target_coords, tuning, seed, control),
    `SQR-DEBIASED-IID` = fit_benchmark_sqr_debiased_iid_v2(train, tau, target_coords, tuning, seed, control),
    `PROFILE-DQR-TRUE-SUPPORT` = fit_benchmark_profile_true_support_v2(
      train, tau, target_coords, tuning, seed,
      active = context$active %||% stop("context$active is required"), control = control
    ),
    `LQMM` = fit_benchmark_lqmm_v2(train, tau, target_coords, tuning, seed,
                                  random_slope = ncol(train$Z) >= 2, control = control),
    `ORACLE-LQMM` = fit_benchmark_lqmm_v2(train, tau, target_coords, tuning, seed,
                                         random_slope = ncol(train$Z) >= 2, control = control),
    `BIAS-ADJ-LQMM` = benchmark_not_implemented_v2(
      method_id, "Faithful two-step/bootstrap implementation from Battagliola et al. is required."
    ),
    `DOUBLE-PEN-QLMM` = benchmark_not_implemented_v2(
      method_id, "Faithful reimplementation of Li, Liu and Luo (2020) is required."
    ),
    `QGEE-SCAD` = benchmark_not_implemented_v2(
      method_id, "Faithful quantile penalised GEE/SCAD implementation from Zu et al. (2023) is required."
    ),
    `QIF-SEE` = benchmark_not_implemented_v2(
      method_id, "Faithful QIF smooth-threshold implementation from Bhattacharya et al. (2026) is required."
    ),
    `PROFILE-DQR-TRUE-NUISANCE` = benchmark_not_implemented_v2(
      method_id, "Oracle score/curvature adapter must use generated random effects."
    ),
    `PROFILE-DQR-POP-H` = benchmark_not_implemented_v2(
      method_id, "Population effective direction must be generated scenario by scenario."
    ),
    `PROFILE-DQR-SPLIT` = benchmark_not_implemented_v2(
      method_id, "Subject-split screening/inference wrapper must be implemented before Module F."
    ),
    `BAYES-MIXED-LASSO` = benchmark_not_implemented_v2(
      method_id, "Legacy mean-model adapter must be relabelled and checked before use."
    ),
    stop("Unknown benchmark method_id: ", method_id)
  )
}
