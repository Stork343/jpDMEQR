# Benchmark adapter layer. Each adapter returns an explicit status and never
# silently substitutes another method.

benchmark_failure_v2 <- function(method_id, stage, message, runtime_sec = NA_real_) {
  list(
    method_id = method_id,
    status = "failed",
    failure_stage = stage,
    failure_message = as.character(message),
    beta_hat = NULL,
    beta_tilde = NULL,
    se = NULL,
    selected = integer(),
    runtime_sec = runtime_sec,
    converged = FALSE,
    warning_messages = character(),
    implementation_version = "sjs-v2"
  )
}

benchmark_not_implemented_v2 <- function(method_id, reason) {
  list(
    method_id = method_id,
    status = "not_implemented",
    failure_stage = "benchmark_not_implemented",
    failure_message = reason,
    beta_hat = NULL,
    beta_tilde = NULL,
    se = NULL,
    selected = integer(),
    runtime_sec = 0,
    converged = FALSE,
    warning_messages = character(),
    implementation_version = "sjs-v2"
  )
}

fit_benchmark_pooled_qr_lasso_v2 <- function(train, tau, target_coords,
                                             tuning, seed, control = list()) {
  if (!requireNamespace("quantreg", quietly = TRUE)) {
    return(benchmark_failure_v2("POOLED-QR-LASSO", "dependency", "quantreg is not installed"))
  }
  set.seed(seed)
  start <- proc.time()[[3]]
  X <- as.matrix(train$X)
  y <- as.numeric(train$y)
  lambda <- tuning$lambda_beta %||% stop("tuning$lambda_beta is required")
  fit <- tryCatch(
    quantreg::rq.fit.lasso(x = X, y = y, tau = tau, lambda = lambda),
    error = function(e) e
  )
  elapsed <- proc.time()[[3]] - start
  if (inherits(fit, "error")) {
    return(benchmark_failure_v2("POOLED-QR-LASSO", "penalised_fit",
                                conditionMessage(fit), elapsed))
  }
  beta <- as.numeric(fit$coefficients)
  names(beta) <- colnames(X) %||% paste0("x", seq_along(beta))
  selected <- which(abs(beta) > (control$selection_tol %||% 1e-8))
  list(
    method_id = "POOLED-QR-LASSO",
    status = "ok",
    failure_stage = "",
    failure_message = "",
    beta_hat = beta,
    beta_tilde = NULL,
    se = NULL,
    selected = selected,
    fit_object = fit,
    prediction_function = function(newX) as.numeric(as.matrix(newX) %*% beta),
    runtime_sec = elapsed,
    converged = TRUE,
    kkt_residual = NA_real_,
    warning_messages = character(),
    implementation_version = "quantreg::rq.fit.lasso"
  )
}

fit_benchmark_profile_dqr_v2 <- function(train, tau, target_coords,
                                         tuning, seed, control = list()) {
  set.seed(seed)
  start <- proc.time()[[3]]
  fit <- tryCatch(
    fit_profile_lasso_v2(
      y = train$y,
      X = train$X,
      Z = train$Z,
      cluster_id = train$cluster_id,
      tau = tau,
      h = tuning$h,
      lambda_beta = tuning$lambda_beta,
      lambda_gamma = tuning$lambda_gamma,
      penalty_factor = control$penalty_factor %||% rep(1, ncol(train$X)),
      control = control$fit_control %||% list()
    ),
    error = function(e) e
  )
  elapsed_fit <- proc.time()[[3]] - start
  if (inherits(fit, "error")) {
    return(benchmark_failure_v2("PROFILE-DQR", "penalised_fit",
                                conditionMessage(fit), elapsed_fit))
  }

  inf <- tryCatch(
    debias_profile_coordinates_v2(
      fit = fit,
      coordinates = target_coords,
      mu_grid = tuning$mu_grid,
      ci_level = control$ci_level %||% 0.95,
      solver_preference = control$solver_preference %||% c("CLARABEL", "ECOS", "SCS")
    ),
    error = function(e) e
  )
  elapsed <- proc.time()[[3]] - start
  if (inherits(inf, "error")) {
    return(benchmark_failure_v2("PROFILE-DQR", "precision_solver",
                                conditionMessage(inf), elapsed))
  }

  selected <- which(abs(fit$beta) > (control$selection_tol %||% 1e-8))
  status <- if (fit$converged && all(inf$table$feasible)) "ok" else "warning"
  warnings <- c(
    if (!fit$converged) "Penalised profile fit did not meet all stopping rules.",
    if (!all(inf$table$feasible)) "One or more Dantzig rows were infeasible."
  )
  list(
    method_id = "PROFILE-DQR",
    status = status,
    failure_stage = "",
    failure_message = "",
    beta_hat = fit$beta,
    beta_tilde = setNames(inf$table$beta_tilde, inf$table$coordinate),
    se = setNames(inf$table$se, inf$table$coordinate),
    ci_lower = setNames(inf$table$ci_lower, inf$table$coordinate),
    ci_upper = setNames(inf$table$ci_upper, inf$table$coordinate),
    selected = selected,
    fit_object = fit,
    inference_object = inf,
    runtime_sec = elapsed,
    converged = fit$converged,
    kkt_residual = fit$kkt_residual,
    profile_identity_error = fit$components$profile_identity_error,
    max_nuisance_gradient = fit$components$max_nuisance_gradient,
    warning_messages = warnings,
    implementation_version = "reference-profile-v2"
  )
}

fit_benchmark_sqr_debiased_iid_v2 <- function(train, tau, target_coords,
                                              tuning, seed, control = list()) {
  iid <- train
  iid$Z <- matrix(0, nrow = length(train$y), ncol = 1,
                  dimnames = list(NULL, "zero_nuisance"))
  iid$cluster_id <- seq_along(train$y)
  ans <- fit_benchmark_profile_dqr_v2(iid, tau, target_coords, tuning, seed, control)
  ans$method_id <- "SQR-DEBIASED-IID"
  ans$implementation_version <- "profile-v2-iid-special-case"
  ans
}

fit_benchmark_profile_true_support_v2 <- function(train, tau, target_coords,
                                                  tuning, seed, active,
                                                  control = list()) {
  target_idx <- if (is.character(target_coords)) match(target_coords, colnames(train$X)) else as.integer(target_coords)
  if (anyNA(target_idx)) stop("Unknown target coordinate in oracle-support adapter.")
  keep <- sort(unique(c(as.integer(active), target_idx)))
  sub <- train
  sub$X <- train$X[, keep, drop = FALSE]
  target_names <- colnames(train$X)[target_idx]
  ans <- fit_benchmark_profile_dqr_v2(
    sub, tau, target_names, tuning, seed,
    control = modifyList(control, list(penalty_factor = rep(0, length(keep))))
  )
  full_beta <- setNames(rep(0, ncol(train$X)), colnames(train$X))
  full_beta[names(ans$beta_hat)] <- ans$beta_hat
  ans$beta_hat <- full_beta
  ans$selected <- keep
  ans$method_id <- "PROFILE-DQR-TRUE-SUPPORT"
  ans$implementation_version <- "profile-v2-oracle-support"
  ans
}
