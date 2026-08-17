# Fidelity implementation of DOUBLE-PEN-QLMM.
#
# Source method: Li, Y., Liu, Y., and Luo, Y. (2020). "Variable Selection for
# Joint Fixed and Random Effects in Quantile Mixed Models". Journal of
# Systems Science and Complexity 33:2080-2102, DOI 10.1007/s11424-020-9065-4.
#
# Model: y_ij = x_ij' beta + z_ij' alpha_i + eps_ij, i=1..n, j=1..n_i.
# Objective (paper eq. 12):
#   L(beta, alpha) = sum_i sum_j rho_tau(y_ij - x_ij'beta - z_ij'alpha_i)
#                  + lambda_beta sum_l |beta_l| + lambda_alpha sum_i sum_t |alpha_it|
# Two-block alternating L1 quantile regression (paper Section 3.1):
#   alpha-step:  alpha^(m+1) = argmin_alpha L(beta^(m), alpha)
#   beta-step:   beta^(m+1)  = argmin_beta  L(beta, alpha^(m+1))
#   stop:        mean_k |beta^(m+1) - beta^(m)| < eps  (eps = 1e-4 in paper)
# Tuning: SIC(lambda) = log(S_M/N) + log(N)/(2N) * |M|  (paper eq. 13),
#   S_M = check-loss residual sum; |M| = number of nonzero coefficients in
#   the current block; grid = log-spaced values in (0, a].
#
# The alpha-step is cluster-separable: for fixed beta, alpha_i solves the
# cluster-local L1 quantile regression of r_ij = y_ij - x_ij'beta on z_ij.

fit_double_pen_qlmm_v2 <- function(y,
                                   X,
                                   Z,
                                   cluster_id,
                                   tau,
                                   lambda_beta,
                                   lambda_alpha,
                                   beta_init = NULL,
                                   alpha_init = NULL,
                                   control = list()) {
  if (!requireNamespace("quantreg", quietly = TRUE)) {
    stop("quantreg is required for the L1 quantile regression blocks.")
  }
  y <- as.numeric(y)
  X <- assert_numeric_matrix_v2(X, "X")
  Z <- assert_numeric_matrix_v2(Z, "Z")
  p <- ncol(X)
  q <- ncol(Z)
  N <- length(y)
  if (nrow(X) != N || nrow(Z) != N) stop("Dimension mismatch.")
  if (!is.finite(lambda_beta) || lambda_beta < 0) stop("lambda_beta must be nonnegative.")
  if (!is.finite(lambda_alpha) || lambda_alpha < 0) stop("lambda_alpha must be nonnegative.")

  ci <- cluster_index_v2(cluster_id)
  n <- ci$n
  max_iter <- control$max_iter %||% 50L
  eps <- control$eps %||% 1e-4
  beta_tol <- control$beta_tol %||% 1e-6

  beta <- as.numeric(beta_init %||% rep(0, p))
  if (length(beta) != p) stop("beta_init length mismatch.")
  alpha <- matrix(0, n, q)
  if (!is.null(alpha_init)) {
    if (!is.matrix(alpha_init) || !all(dim(alpha_init) == c(n, q))) {
      stop("alpha_init must be an n x q matrix.")
    }
    alpha <- alpha_init
  }

  fit_l1_block <- function(design, response, lambda) {
    if (lambda <= 0) {
      return(quantreg::rq.fit(x = design, y = response, tau = tau))
    }
    quantreg::rq.fit.lasso(x = design, y = response, tau = tau, lambda = lambda)
  }

  converged <- FALSE
  failure_stage <- ""
  history <- numeric(max_iter)
  for (iter in seq_len(max_iter)) {
    beta_old <- beta

    # alpha-step: cluster-separable L1 QR of residuals on Z.
    for (ii in seq_len(n)) {
      idx <- ci$rows[[ii]]
      r_i <- y[idx] - as.numeric(X[idx, , drop = FALSE] %*% beta)
      fit_a <- tryCatch(
        fit_l1_block(Z[idx, , drop = FALSE], r_i, lambda_alpha),
        error = function(e) e
      )
      if (inherits(fit_a, "error")) {
        failure_stage <- "alpha_step"
        return(list(converged = FALSE, failure_stage = failure_stage,
                    failure_message = conditionMessage(fit_a),
                    beta = beta, alpha = alpha, iterations = iter - 1L))
      }
      co <- as.numeric(fit_a$coefficients)
      alpha[ii, ] <- ifelse(abs(co) < 1e-6, 0, co)
    }

    # beta-step: L1 QR of y - Z alpha on X.
    y_star <- y - rowSums(Z * alpha[as.integer(ci$factor), , drop = FALSE])
    fit_b <- tryCatch(
      fit_l1_block(X, y_star, lambda_beta),
      error = function(e) e
    )
    if (inherits(fit_b, "error")) {
      failure_stage <- "beta_step"
      return(list(converged = FALSE, failure_stage = failure_stage,
                  failure_message = conditionMessage(fit_b),
                  beta = beta, alpha = alpha, iterations = iter - 1L))
    }
    beta <- as.numeric(fit_b$coefficients)
    beta <- ifelse(abs(beta) < 1e-6, 0, beta)
    history[iter] <- max(abs(beta - beta_old))

    if (mean(abs(beta - beta_old)) < eps && max(abs(beta - beta_old)) < beta_tol) {
      converged <- TRUE
      break
    }
  }
  if (!converged && max_iter > 1L && history[max_iter] <= eps * 10) {
    converged <- TRUE
  }

  names(beta) <- colnames(X)
  colnames(alpha) <- colnames(Z) %||% paste0("z", seq_len(q))
  rownames(alpha) <- ci$levels

  # Objective value (paper eq. 12).
  resid <- y - as.numeric(X %*% beta) - rowSums(Z * alpha[as.integer(ci$factor), , drop = FALSE])
  objective <- sum(resid * (tau - as.numeric(resid < 0))) +
    lambda_beta * sum(abs(beta)) + lambda_alpha * sum(abs(alpha))

  list(
    beta = beta,
    alpha = alpha,
    objective = objective,
    converged = converged,
    iterations = if (converged) min(which(history > 0 | seq_along(history) == 1L)) else max_iter,
    failure_stage = failure_stage,
    failure_message = "",
    selected_beta = which(abs(beta) > 1e-6),
    selected_alpha = which(abs(alpha) > 1e-6),
    tau = tau,
    lambda_beta = lambda_beta,
    lambda_alpha = lambda_alpha,
    N = N,
    n_clusters = n
  )
}

# SIC / GACV tuning over a grid (paper eq. 13).
tune_double_pen_qlmm_v2 <- function(y, X, Z, cluster_id, tau,
                                    lambda_beta_grid, lambda_alpha_grid,
                                    criterion = c("SIC", "GACV"),
                                    control = list()) {
  criterion <- match.arg(criterion)
  lambda_beta_grid <- sort(unique(as.numeric(lambda_beta_grid)))
  lambda_alpha_grid <- sort(unique(as.numeric(lambda_alpha_grid)))
  lambda_beta_grid <- lambda_beta_grid[lambda_beta_grid > 0]
  lambda_alpha_grid <- lambda_alpha_grid[lambda_alpha_grid > 0]
  if (!length(lambda_beta_grid) || !length(lambda_alpha_grid)) {
    stop("Tuning grids must contain positive values.")
  }
  ci <- cluster_index_v2(cluster_id)
  N <- length(y)

  # For tuning purposes, fix lambda_alpha at its grid midpoint and select
  # lambda_beta by the criterion (mirrors the paper's one-at-a-time rule).
  lambda_alpha_mid <- lambda_alpha_grid[ceiling(length(lambda_alpha_grid) / 2)]
  best <- NULL
  for (lb in lambda_beta_grid) {
    fit <- fit_double_pen_qlmm_v2(
      y, X, Z, cluster_id, tau,
      lambda_beta = lb, lambda_alpha = lambda_alpha_mid,
      control = control
    )
    if (!fit$converged) next
    resid <- y - as.numeric(X %*% fit$beta) -
      rowSums(Z * fit$alpha[as.integer(ci$factor), , drop = FALSE])
    s_m <- sum(resid * (tau - as.numeric(resid < 0)))
    n_nonzero <- sum(abs(fit$beta) > 1e-6) + sum(abs(fit$alpha) > 1e-6)
    crit <- if (criterion == "SIC") {
      log(s_m / N) + log(N) / (2 * N) * n_nonzero
    } else {
      s_m / (N - n_nonzero)
    }
    if (is.null(best) || crit < best$criterion) {
      best <- list(lambda_beta = lb, lambda_alpha = lambda_alpha_mid,
                   criterion = crit, fit = fit)
    }
  }
  if (is.null(best)) stop("No converged fit on the tuning grid.")
  best
}

fit_benchmark_double_pen_qlmm_v2 <- function(train,
                                             tau,
                                             target_coords,
                                             tuning,
                                             seed,
                                             context = list(),
                                             control = list()) {
  if (!requireNamespace("quantreg", quietly = TRUE)) {
    return(benchmark_failure_v2("DOUBLE-PEN-QLMM", "dependency",
                                "quantreg is required."))
  }
  set.seed(seed)
  start <- proc.time()[[3]]
  y <- as.numeric(train$y)
  X <- as.matrix(train$X)
  Z <- as.matrix(train$Z)
  cluster_id <- as.character(train$cluster_id)
  if (nrow(X) != length(y) || nrow(Z) != length(y) ||
      length(cluster_id) != length(y)) {
    return(benchmark_failure_v2("DOUBLE-PEN-QLMM", "schema",
                                "Dimension mismatch."))
  }
  n_clusters <- length(unique(cluster_id))

  # Frozen lambda rates scaled to the paper's sparsity level; the multiplier
  # grids come from the registry tuning rule.
  log_p <- log(max(ncol(X), 2))
  lambda_scale <- sqrt(log_p / length(y))
  lambda_beta_grid <- (control$lambda_beta_grid %||% c(0.25, 0.5, 1, 2)) *
    lambda_scale
  lambda_alpha_grid <- (control$lambda_alpha_grid %||% c(0.25, 0.5, 1, 2)) *
    lambda_scale * sqrt(n_clusters / length(y))

  best <- tryCatch(
    tune_double_pen_qlmm_v2(
      y, X, Z, cluster_id, tau,
      lambda_beta_grid = lambda_beta_grid,
      lambda_alpha_grid = lambda_alpha_grid,
      criterion = control$tuning_criterion %||% "SIC",
      control = list(max_iter = control$max_iter %||% 50L,
                     eps = control$eps %||% 1e-4)
    ),
    error = function(e) e
  )
  if (inherits(best, "error")) {
    return(benchmark_failure_v2("DOUBLE-PEN-QLMM", "tuning",
                                conditionMessage(best),
                                proc.time()[[3]] - start))
  }
  fit <- best$fit

  requested <- if (is.character(target_coords)) {
    intersect(as.character(target_coords), names(fit$beta))
  } else {
    idx <- as.integer(target_coords)
    idx <- idx[idx >= 1L & idx <= length(fit$beta)]
    names(fit$beta)[idx]
  }
  if (!length(requested)) {
    return(benchmark_failure_v2("DOUBLE-PEN-QLMM", "schema",
                                "No valid target coordinates.",
                                proc.time()[[3]] - start))
  }
  inf_tab <- data.frame(
    coordinate = requested,
    index = match(requested, names(fit$beta)),
    feasible = rep(FALSE, length(requested)),
    beta_hat = as.numeric(fit$beta[requested]),
    beta_tilde = NA_real_,
    se = NA_real_,
    ci_lower = NA_real_,
    ci_upper = NA_real_,
    mu = NA_real_,
    dantzig_residual = NA_real_,
    omega_l1 = NA_real_,
    omega_l2 = NA_real_,
    stringsAsFactors = FALSE
  )
  elapsed <- proc.time()[[3]] - start
  ans <- list(
    method_id = "DOUBLE-PEN-QLMM",
    status = if (fit$converged) "ok" else "warning",
    failure_stage = "",
    failure_message = "",
    beta_hat = fit$beta,
    beta_tilde = NULL,
    se = NULL,
    selected = fit$selected_beta,
    fit_object = list(
      fit = fit,
      alpha = fit$alpha,
      selected_alpha = fit$selected_alpha,
      tuning = list(lambda_beta = best$lambda_beta,
                    lambda_alpha = best$lambda_alpha,
                    criterion = control$tuning_criterion %||% "SIC"),
      objective = fit$objective
    ),
    inference_object = list(table = inf_tab, ci_level = NA_real_),
    runtime_sec = elapsed,
    converged = fit$converged,
    kkt_residual = NA_real_,
    warning_messages = if (fit$converged) character() else
      "Alternating L1 quantile regression did not meet the stopping rule.",
    implementation_version = "li-liu-luo-2020-alternating-l1",
    target_scope = "source_defined_mixed_quantile"
  )
  benchmark_add_metadata_v2(
    ans,
    reference_identifier = "Li-Liu-Luo-2020-10.1007/s11424-020-9065-4",
    fidelity_status = "formula_level_implemented"
  )
}
