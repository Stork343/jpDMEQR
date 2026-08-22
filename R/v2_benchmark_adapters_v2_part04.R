# Formula-level fidelity implementation of SQR-DEBIASED-IID.
#
# Source method: Yan, Wang and Zhang (2023), "Confidence Intervals and
# Hypothesis Testing for High-dimensional Quantile Regression: Convolution
# Smoothing and Debiasing", JMLR 24(245):1-49, paper 22-1217.
#
# This adapter exposes the iid convolution-smoothed quantile objective,
# iid score/Hessian, row-wise CLIME precision direction, one-step debiased
# update and iid variance formula directly. It deliberately does NOT reuse
# the profiled nuisance machinery of PROFILE-DQR: there is no nuisance
# ridge, no cluster sandwich and no profiling in this adapter.
#
# Formula inventory (JMLR 22-1217):
#   l_{h,tau}(u)   = (rho_tau * K_h)(u);  Gaussian kernel: u*(tau - Phi(-u/h)) + h*phi(u/h)
#   score(beta)    = n^{-1} sum_i [ Kbar((x_i'beta - y_i)/h) - tau ] x_i
#   Hessian(beta)  = (n h)^{-1} sum_i K((x_i'beta - y_i)/h) x_i x_i'
#   l1-SQR         = argmin n^{-1} sum_i l_{h,tau}(y_i - x_i'beta) + lambda ||beta||_1
#   W_hat          = argmin ||W||_1 s.t. || W Hessian - I ||_inf <= gamma   (row-wise)
#   debiased       = beta_hat + n^{-1} W_hat sum_i [ tau - Kbar((x_i'beta_hat - y_i)/h) ] x_i
#   variance       = tau(1-tau) * w' Sigma_hat w / N,  Sigma_hat = N^{-1} sum x_i x_i'
#   bandwidth      = max{ sqrt(tau(1-tau)) (log p)^{1/2} N^{-3/10}, 0.05 }
#   lambda rate    ~ sqrt( tau(1-tau) log p / N )
#
# The observation-level sample size N (not the cluster count) is the
# stochastic index of the iid estimator and variance, matching the source.

# ---------------------------------------------------------------------------
# Kernel and smoothed loss (Gaussian kernel, JMLR eq. 6)
# ---------------------------------------------------------------------------

sqr_iid_smooth_loss_v2 <- function(u, tau, h) {
  u * (tau - stats::pnorm(-u / h)) + h * stats::dnorm(u / h)
}

sqr_iid_score_fn_v2 <- function(y, X, beta, tau, h) {
  r <- y - as.numeric(X %*% beta)
  # JMLR eq.: score = n^{-1} sum_i [ Kbar((x_i'beta - y_i)/h) - tau ] x_i.
  # Here (x'beta - y)/h = -r/h, and Phi(-r/h) = 1 - Phi(r/h).
  as.numeric(crossprod(X, 1 - stats::pnorm(r / h) - tau)) / length(y)
}

sqr_iid_hessian_fn_v2 <- function(y, X, beta, tau, h) {
  r <- y - as.numeric(X %*% beta)
  w <- stats::dnorm(r / h) / h
  crossprod(X, X * w) / length(y)
}

# ---------------------------------------------------------------------------
# l1-penalised smoothed quantile regression by proximal gradient
# (JMLR Algorithm 1, step 1). Deterministic, dimension-validated.
# ---------------------------------------------------------------------------

fit_sqr_iid_lasso_v2 <- function(y,
                                 X,
                                 tau,
                                 h,
                                 lambda,
                                 penalty_factor = NULL,
                                 beta_init = NULL,
                                 control = list()) {
  y <- as.numeric(y)
  X <- assert_numeric_matrix_v2(X, "X")
  tau <- as.numeric(tau)
  h <- as.numeric(h)
  if (tau <= 0 || tau >= 1) stop("tau must lie in (0,1).")
  if (!is.finite(h) || h <= 0) stop("h must be positive.")
  if (!is.finite(lambda) || lambda < 0) stop("lambda must be nonnegative.")
  p <- ncol(X)
  N <- length(y)
  if (nrow(X) != N) stop("y and X dimension mismatch.")

  penalty_factor <- as.numeric(penalty_factor %||% rep(1, p))
  if (length(penalty_factor) != p || any(!is.finite(penalty_factor)) ||
      any(penalty_factor < 0)) stop("Invalid penalty_factor.")

  beta <- as.numeric(beta_init %||% rep(0, p))
  if (length(beta) != p) stop("beta_init length mismatch.")

  max_iter <- control$max_iter %||% 500L
  beta_tol <- control$beta_tol %||% 1e-7
  kkt_tol <- control$kkt_tol %||% 1e-5
  backtrack_factor <- control$backtrack_factor %||% 0.5
  max_backtrack <- control$max_backtrack %||% 40L
  verbose <- isTRUE(control$verbose)

  objective <- function(b) {
    mean(sqr_iid_smooth_loss_v2(y - as.numeric(X %*% b), tau, h)) +
      lambda * sum(penalty_factor * abs(b))
  }
  smooth_value <- function(b) {
    mean(sqr_iid_smooth_loss_v2(y - as.numeric(X %*% b), tau, h))
  }
  score0 <- sqr_iid_score_fn_v2(y, X, beta, tau, h)
  value <- objective(beta)
  converged <- FALSE
  failure_stage <- ""

  for (iter in seq_len(max_iter)) {
    # Local Lipschitz estimate from the largest eigenvalue of the Hessian.
    H <- sqr_iid_hessian_fn_v2(y, X, beta, tau, h)
    eigmax <- largest_eigenvalue_power_v2(H, maxit = 100L, tol = 1e-8)
    step <- 1 / max(eigmax, 1e-6)
    score <- sqr_iid_score_fn_v2(y, X, beta, tau, h)
    accepted <- FALSE
    candidate <- NULL
    f_beta <- smooth_value(beta)

    for (bt in 0:max_backtrack) {
      beta_new <- soft_threshold_weighted_v2(
        beta - step * score,
        step * lambda * penalty_factor
      )
      delta <- beta_new - beta
      cand_value <- objective(beta_new)
      # Correct proximal-gradient majorant: the smooth part is majorised at
      # beta, while the l1 penalty is evaluated at the candidate beta_new.
      upper <- f_beta + sum(score * delta) + sum(delta^2) / (2 * step) +
        lambda * sum(penalty_factor * abs(beta_new)) + 1e-12
      if (cand_value <= upper) {
        accepted <- TRUE
        candidate <- beta_new
        break
      }
      step <- step * backtrack_factor
      if (step < 1e-14) break
    }

    if (!accepted) {
      failure_stage <- "sqr_backtracking"
      break
    }
    beta <- candidate
    value <- objective(beta)
    # Converge only when BOTH the parameter change and the KKT residual are
    # small; a small first step from the zero start must not be mistaken for
    # convergence.
    score_check <- sqr_iid_score_fn_v2(y, X, beta, tau, h)
    kkt_check <- max(abs(kkt_residual_profile_v2(
      beta, score_check, lambda, penalty_factor
    )))
    if (max(abs(beta_new - beta)) < beta_tol * max(1, max(abs(beta))) &&
        kkt_check <= kkt_tol) {
      converged <- TRUE
      break
    }
    if (verbose && iter %% 50L == 0L) {
      message(sprintf("sqr-iid iter %d value %.8f", iter, value))
    }
  }

  score_final <- sqr_iid_score_fn_v2(y, X, beta, tau, h)
  kkt <- kkt_residual_profile_v2(
    beta, score_final, lambda, penalty_factor
  )
  if (!converged && max(abs(kkt)) <= kkt_tol) converged <- TRUE
  names(beta) <- colnames(X)

  list(
    beta = beta,
    objective = value,
    score = score_final,
    hessian = sqr_iid_hessian_fn_v2(y, X, beta, tau, h),
    kkt_residual = max(abs(kkt)),
    converged = converged,
    iterations = iter,
    failure_stage = failure_stage,
    h = h,
    lambda = lambda,
    tau = tau,
    N = N
  )
}

# ---------------------------------------------------------------------------
# Row-wise CLIME precision direction (JMLR eq. 11/15) via CVXR.
# ---------------------------------------------------------------------------

solve_clime_row_iid_v2 <- function(H,
                                   coordinate,
                                   gamma_grid,
                                   solver_preference = c("CLARABEL", "ECOS", "SCS"),
                                   solver_opts = list()) {
  if (!requireNamespace("CVXR", quietly = TRUE)) {
    stop("Package 'CVXR' is required for the CLIME row program.")
  }
  H <- assert_numeric_matrix_v2(H, "H")
  H <- (H + t(H)) / 2
  p <- nrow(H)
  if (ncol(H) != p) stop("H must be square.")

  k <- if (is.character(coordinate)) {
    if (is.null(colnames(H)) || !coordinate %in% colnames(H)) {
      stop("Unknown coordinate name: ", coordinate)
    }
    match(coordinate, colnames(H))
  } else {
    as.integer(coordinate)
  }
  if (length(k) != 1L || k < 1L || k > p) stop("Invalid coordinate index.")

  gamma_grid <- sort(unique(as.numeric(gamma_grid)))
  if (!length(gamma_grid) || any(!is.finite(gamma_grid)) || any(gamma_grid <= 0)) {
    stop("gamma_grid must contain positive finite values.")
  }
  installed <- CVXR::installed_solvers()
  solver <- solver_preference[solver_preference %in% installed][1]
  if (is.na(solver)) stop("No requested CVXR solver is installed.")

  cvxr_solve <- if ("psolve" %in% getNamespaceExports("CVXR")) {
    getExportedValue("CVXR", "psolve")
  } else {
    getExportedValue("CVXR", "solve")
  }

  e <- numeric(p); e[k] <- 1
  attempts <- vector("list", length(gamma_grid))
  for (ii in seq_along(gamma_grid)) {
    gamma <- gamma_grid[ii]
    b_var <- CVXR::Variable(p)
    constraints <- list(
      H %*% b_var - e <= gamma,
      H %*% b_var - e >= -gamma
    )
    problem <- CVXR::Problem(
      CVXR::Minimize(CVXR::p_norm(b_var, 1)),
      constraints
    )
    solved <- tryCatch(
      do.call(cvxr_solve, c(list(problem = problem, solver = solver), solver_opts)),
      error = function(e) e
    )
    if (inherits(solved, "error")) {
      attempts[[ii]] <- data.frame(gamma = gamma, status = "solver_error",
                                   residual = NA_real_, message = conditionMessage(solved))
      next
    }
    if (is.numeric(solved) && length(solved) == 1L &&
        "value" %in% getNamespaceExports("CVXR") &&
        "status" %in% getNamespaceExports("CVXR")) {
      # CVXR >= 1.9 contract: solved is the objective value.
      status <- tryCatch(CVXR::status(problem), error = function(e) "unknown")
      omega <- tryCatch(as.numeric(CVXR::value(b_var)), error = function(e) NULL)
    } else {
      status <- solved$status %||% "unknown"
      omega <- tryCatch(as.numeric(solved$getValue(b_var)), error = function(e) NULL)
    }
    if (is.null(omega) || length(omega) != p || any(!is.finite(omega))) {
      attempts[[ii]] <- data.frame(gamma = gamma, status = status,
                                   residual = NA_real_, message = "No finite solution")
      next
    }
    residual <- max(abs(as.numeric(H %*% omega - e)))
    attempts[[ii]] <- data.frame(gamma = gamma, status = status,
                                 residual = residual, message = "")
    if (residual <= gamma + 1e-6 && status %in% c("optimal", "optimal_inaccurate")) {
      return(list(
        omega = omega,
        coordinate = k,
        coordinate_name = colnames(H)[k] %||% as.character(k),
        gamma = gamma,
        residual = residual,
        solver = solver,
        solver_status = status,
        l1_norm = sum(abs(omega)),
        l2_norm = sqrt(sum(omega^2)),
        attempts = do.call(rbind, attempts[seq_len(ii)]),
        feasible = TRUE
      ))
    }
  }
  list(
    omega = rep(NA_real_, p),
    coordinate = k,
    coordinate_name = colnames(H)[k] %||% as.character(k),
    gamma = NA_real_,
    residual = NA_real_,
    solver = solver,
    solver_status = "infeasible_on_grid",
    l1_norm = NA_real_,
    l2_norm = NA_real_,
    attempts = do.call(rbind, attempts),
    feasible = FALSE
  )
}

# ---------------------------------------------------------------------------
# Benchmark adapter with the required fit_benchmark_<id>_v2 signature.
# ---------------------------------------------------------------------------

fit_benchmark_sqr_debiased_iid_v2 <- function(train,
                                              tau,
                                              target_coords,
                                              tuning,
                                              seed,
                                              control = list()) {
  set.seed(seed)
  start <- proc.time()[[3]]
  y <- as.numeric(train$y)
  X <- as.matrix(train$X)
  N <- length(y)
  p <- ncol(X)
  if (nrow(X) != N) {
    return(benchmark_failure_v2("SQR-DEBIASED-IID", "schema",
                                "y and X dimension mismatch."))
  }

  # Source-method bandwidth and penalty rates (JMLR 22-1217 Section 4).
  h_iid <- max(sqrt(tau * (1 - tau)) * sqrt(log(max(p, 2))) * N^(-3 / 10), 0.05)
  lambda_mult <- if (!is.null(tuning$lambda_beta_multiplier)) {
    tuning$lambda_beta_multiplier
  } else if (!is.null(tuning$lambda_beta) && is.finite(tuning$lambda_beta) &&
             tuning$lambda_beta > 0) {
    1
  } else {
    1
  }
  lambda <- lambda_mult * sqrt(tau * (1 - tau) * log(max(p, 2)) / N)
  gamma_grid <- if (!is.null(tuning$mu_grid) && length(tuning$mu_grid)) {
    as.numeric(tuning$mu_grid)
  } else {
    sqrt(log(max(p, 2)) / (N * h_iid)) + h_iid^2
  }

  fit <- tryCatch(
    fit_sqr_iid_lasso_v2(
      y, X, tau, h = h_iid, lambda = lambda,
      penalty_factor = control$penalty_factor %||% rep(1, p),
      control = control$fit_control %||% list()
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(benchmark_failure_v2("SQR-DEBIASED-IID", "penalised_fit",
                                conditionMessage(fit), proc.time()[[3]] - start))
  }

  coords <- if (is.character(target_coords)) {
    match(target_coords, colnames(X))
  } else {
    as.integer(target_coords)
  }
  coords <- coords[!is.na(coords)]
  if (!length(coords)) {
    return(benchmark_failure_v2("SQR-DEBIASED-IID", "schema",
                                "No valid target coordinates."))
  }

  alpha <- 1 - (control$ci_level %||% 0.95)
  zcrit <- stats::qnorm(1 - alpha / 2)
  Sigma_hat <- crossprod(X) / N
  rows <- vector("list", length(coords))
  directions <- vector("list", length(coords))
  names(directions) <- colnames(X)[coords]

  for (ii in seq_along(coords)) {
    k <- coords[ii]
    direction <- tryCatch(
      solve_clime_row_iid_v2(
        fit$hessian, coordinate = k, gamma_grid = gamma_grid,
        solver_preference = control$solver_preference %||% c("CLARABEL", "ECOS", "SCS"),
        solver_opts = control$solver_opts %||% list()
      ),
      error = function(e) e
    )
    if (inherits(direction, "error")) {
      return(benchmark_failure_v2("SQR-DEBIASED-IID", "precision_solver",
                                  conditionMessage(direction), proc.time()[[3]] - start))
    }
    directions[[ii]] <- direction
    if (!direction$feasible) {
      rows[[ii]] <- data.frame(
        coordinate = colnames(X)[k], index = k, feasible = FALSE,
        beta_hat = fit$beta[k], beta_tilde = NA_real_, se = NA_real_,
        ci_lower = NA_real_, ci_upper = NA_real_, mu = NA_real_,
        dantzig_residual = NA_real_, omega_l1 = NA_real_, omega_l2 = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }
    omega <- direction$omega
    # One-step debiasing (JMLR eq. 10):
    #   beta_tilde = beta_hat + n^{-1} W sum_i [ tau - Kbar((x_i'beta_hat - y_i)/h) ] x_i
    #             = beta_hat - W * score,  since score = n^{-1} sum [Kbar(.) - tau] x.
    beta_tilde <- as.numeric(fit$beta[k] - sum(omega * fit$score))
    var_k <- tau * (1 - tau) * drop(crossprod(omega, Sigma_hat %*% omega)) / N
    se <- sqrt(max(var_k, 0))
    rows[[ii]] <- data.frame(
      coordinate = colnames(X)[k],
      index = k,
      feasible = TRUE,
      beta_hat = as.numeric(fit$beta[k]),
      beta_tilde = beta_tilde,
      se = se,
      ci_lower = beta_tilde - zcrit * se,
      ci_upper = beta_tilde + zcrit * se,
      mu = direction$gamma,
      dantzig_residual = direction$residual,
      omega_l1 = direction$l1_norm,
      omega_l2 = direction$l2_norm,
      stringsAsFactors = FALSE
    )
  }

  tab <- do.call(rbind, rows)
  elapsed <- proc.time()[[3]] - start
  ans <- list(
    method_id = "SQR-DEBIASED-IID",
    status = if (fit$converged) "ok" else "warning",
    failure_stage = "",
    failure_message = "",
    beta_hat = fit$beta,
    beta_tilde = setNames(tab$beta_tilde, tab$coordinate),
    se = setNames(tab$se, tab$coordinate),
    ci_lower = setNames(tab$ci_lower, tab$coordinate),
    ci_upper = setNames(tab$ci_upper, tab$coordinate),
    selected = which(abs(fit$beta) > (control$selection_tol %||% 1e-8)),
    fit_object = fit,
    inference_object = list(
      table = tab,
      directions = directions,
      ci_level = control$ci_level %||% 0.95,
      iid = TRUE,
      variance_formula = "tau(1-tau) * omega' Sigma_hat omega / N",
      bandwidth = h_iid,
      lambda = lambda,
      gamma_grid = gamma_grid
    ),
    runtime_sec = elapsed,
    converged = fit$converged,
    kkt_residual = fit$kkt_residual,
    warning_messages = if (fit$converged) character() else
      "Smoothed QR fit did not meet all stopping rules.",
    implementation_version = "sqr-iid-fidelity-jmlr22-1217",
    target_scope = "iid_smoothed_quantile",
    variance_scope = "naive_iid"
  )
  benchmark_add_metadata_v2(
    ans,
    reference_identifier = "Yan-Wang-Zhang-2023-JMLR-22-1217",
    fidelity_status = "formula_level_implemented"
  )
}
