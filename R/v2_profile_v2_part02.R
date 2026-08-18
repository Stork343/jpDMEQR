profile_components_v2 <- function(y,
                                  X,
                                  Z,
                                  cluster_id,
                                  beta,
                                  tau,
                                  h,
                                  lambda_gamma = 1,
                                  gamma_init = NULL,
                                  need_hessian = TRUE,
                                  nuisance_control = list()) {
  validate_tau_h_v2(tau, h)
  y <- as.numeric(y)
  X <- assert_numeric_matrix_v2(X, "X")
  Z <- assert_numeric_matrix_v2(Z, "Z")
  beta <- as.numeric(beta)
  N <- length(y)
  p <- ncol(X)
  q <- ncol(Z)
  if (nrow(X) != N || nrow(Z) != N || length(cluster_id) != N) {
    stop("y, X, Z and cluster_id dimensions do not match.")
  }
  if (length(beta) != p) stop("beta has wrong length.")

  Lambda <- make_lambda_matrix_v2(lambda_gamma, q)
  ci <- cluster_index_v2(cluster_id)
  n <- ci$n

  cluster_scores <- matrix(NA_real_, nrow = n, ncol = p,
                           dimnames = list(ci$levels, colnames(X)))
  H_sum <- if (need_hessian) matrix(0, p, p) else NULL
  H_identity_sum <- if (need_hessian) matrix(0, p, p) else NULL
  gamma_list <- setNames(vector("list", n), ci$levels)
  residual_list <- setNames(vector("list", n), ci$levels)
  nuisance_diagnostics <- vector("list", n)
  objective_values <- numeric(n)
  identity_errors <- numeric(n)

  # Per-cluster worker: pure function of (y_i, X_i, Z_i, beta, ...) with no
  # shared state, so parallel execution is numerically identical to the
  # sequential loop. Each cluster's nuisance solve is deterministic (BFGS +
  # damped Newton; no global RNG dependence beyond the supplied seed).
  cluster_worker <- function(ii) {
    idx <- ci$rows[[ii]]
    cname <- ci$levels[ii]
    y_i <- y[idx]
    X_i <- X[idx, , drop = FALSE]
    Z_i <- Z[idx, , drop = FALSE]
    m_i <- length(idx)
    g0 <- resolve_gamma_init_v2(gamma_init, cname, ii, q)

    ng <- profile_cluster_gamma_v2(
      y_i = y_i, X_i = X_i, Z_i = Z_i, beta = beta,
      tau = tau, h = h, Lambda = Lambda,
      gamma_init = g0, control = nuisance_control
    )
    gamma_i <- ng$gamma
    r_i <- y_i - as.numeric(X_i %*% beta) - as.numeric(Z_i %*% gamma_i)
    psi_i <- psi_smooth_v2(r_i, tau, h)

    diag_row <- data.frame(
      cluster = cname,
      m_i = m_i,
      converged = ng$converged,
      gradient_max = ng$gradient_max,
      optim_convergence = ng$optim_convergence,
      function_evaluations = ng$function_evaluations,
      gradient_evaluations = ng$gradient_evaluations,
      stringsAsFactors = FALSE
    )
    objective_i <- mean(rho_smooth_v2(r_i, tau, h)) +
      0.5 * drop(crossprod(gamma_i, Lambda %*% gamma_i))
    score_i <- as.numeric(-crossprod(X_i, psi_i) / m_i)

    H_i <- NULL
    H_identity_i <- NULL
    identity_err_i <- NA_real_
    if (need_hessian) {
      w <- phi_smooth_v2(r_i, tau, h) / m_i
      WX <- X_i * w
      WZ <- Z_i * w
      H_xx <- crossprod(X_i, WX)
      H_xz <- crossprod(X_i, WZ)
      B_i <- crossprod(Z_i, WZ) + Lambda
      ZWX <- crossprod(Z_i, WX)
      A_i <- solve_linear_pd_v2(B_i, ZWX, name = paste0("B_i for cluster ", cname))
      H_i <- H_xx - H_xz %*% A_i

      X_tilde <- X_i - Z_i %*% A_i
      H_identity_i <- crossprod(X_tilde, X_tilde * w) +
        t(A_i) %*% Lambda %*% A_i
      identity_err_i <- max(abs(H_i - H_identity_i))
    }
    list(
      ii = ii,
      gamma = gamma_i,
      residual = r_i,
      diag = diag_row,
      objective = objective_i,
      score = score_i,
      H = H_i,
      H_identity = H_identity_i,
      identity_error = identity_err_i
    )
  }

  cluster_cores <- getOption("jpDMEQR.cluster_cores", 1L)
  out <- if (cluster_cores > 1L && .Platform$OS.type != "windows" && n > 1L) {
    parallel::mclapply(seq_len(n), cluster_worker,
                       mc.cores = min(cluster_cores, n),
                       mc.preschedule = FALSE)
  } else {
    lapply(seq_len(n), cluster_worker)
  }
  for (item in out) {
    ii <- item$ii
    gamma_list[[ii]] <- item$gamma
    residual_list[[ii]] <- item$residual
    nuisance_diagnostics[[ii]] <- item$diag
    objective_values[ii] <- item$objective
    cluster_scores[ii, ] <- item$score
    if (need_hessian) {
      H_sum <- H_sum + item$H
      H_identity_sum <- H_identity_sum + item$H_identity
      identity_errors[ii] <- item$identity_error
    }
  }

  nuisance_df <- do.call(rbind, nuisance_diagnostics)
  gamma_mat <- do.call(rbind, gamma_list)
  rownames(gamma_mat) <- ci$levels
  colnames(gamma_mat) <- colnames(Z) %||% paste0("z", seq_len(q))

  score <- colMeans(cluster_scores)
  hessian <- if (need_hessian) (H_sum / n + t(H_sum / n)) / 2 else NULL
  hessian_identity <- if (need_hessian) (H_identity_sum / n + t(H_identity_sum / n)) / 2 else NULL

  list(
    objective = mean(objective_values),
    score = as.numeric(score),
    hessian = hessian,
    hessian_identity = hessian_identity,
    profile_identity_error = if (need_hessian) max(abs(hessian - hessian_identity)) else NA_real_,
    cluster_identity_errors = identity_errors,
    cluster_scores = cluster_scores,
    gamma = gamma_mat,
    gamma_list = gamma_list,
    residuals = residual_list,
    nuisance_diagnostics = nuisance_df,
    nuisance_converged = all(nuisance_df$converged),
    max_nuisance_gradient = max(nuisance_df$gradient_max),
    Lambda = Lambda,
    n_clusters = n,
    total_rows = N,
    p = p,
    q = q,
    tau = tau,
    h = h
  )
}

profile_objective_v2 <- function(y, X, Z, cluster_id, beta, tau, h,
                                 lambda_gamma = 1, gamma_init = NULL,
                                 nuisance_control = list()) {
  profile_components_v2(
    y = y, X = X, Z = Z, cluster_id = cluster_id,
    beta = beta, tau = tau, h = h, lambda_gamma = lambda_gamma,
    gamma_init = gamma_init, need_hessian = FALSE,
    nuisance_control = nuisance_control
  )$objective
}

# -----------------------------------------------------------------------------
# Penalised profile estimator
# -----------------------------------------------------------------------------

soft_threshold_weighted_v2 <- function(x, threshold) {
  sign(x) * pmax(abs(x) - threshold, 0)
}

largest_eigenvalue_power_v2 <- function(A, maxit = 200L, tol = 1e-8) {
  A <- (A + t(A)) / 2
  p <- nrow(A)
  if (p == 1L) return(max(A[1, 1], 0))
  v <- rep(1 / sqrt(p), p)
  value_old <- 0
  for (iter in seq_len(maxit)) {
    Av <- as.numeric(A %*% v)
    nv <- sqrt(sum(Av^2))
    if (!is.finite(nv) || nv == 0) return(0)
    v <- Av / nv
    value <- drop(crossprod(v, A %*% v))
    if (abs(value - value_old) <= tol * (1 + abs(value))) break
    value_old <- value
  }
  max(as.numeric(value), 0)
}

kkt_residual_profile_v2 <- function(beta, score, lambda_beta, penalty_factor,
                                    active_tol = 1e-10) {
  beta <- as.numeric(beta)
  score <- as.numeric(score)
  pf <- as.numeric(penalty_factor)
  out <- numeric(length(beta))
  unpen <- pf == 0
  out[unpen] <- abs(score[unpen])

  pen <- !unpen
  active <- pen & abs(beta) > active_tol
  inactive <- pen & !active
  out[active] <- abs(score[active] + lambda_beta * pf[active] * sign(beta[active]))
  out[inactive] <- pmax(0, abs(score[inactive]) - lambda_beta * pf[inactive])
  max(out)
}
