# -----------------------------------------------------------------------------
# Unsmoothed quantile score (primary Wald meat)
# -----------------------------------------------------------------------------

# psi_tau(u) = tau - I(u < 0), with the deterministic convention I(u < 0) so
# that psi(0) = tau. Used for the primary cluster meat of the unsmoothed
# regularised profile target (theory decision: PILOT_GATE_THEORY_DECISION.md).
psi_quantile_v2 <- function(u, tau) {
  tau - as.numeric(u < 0)
}

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
  unsmoothed_cluster_scores <- cluster_scores
  H_sum <- if (need_hessian) matrix(0, p, p) else NULL
  H_identity_sum <- if (need_hessian) matrix(0, p, p) else NULL
  hessian_fold_sums <- if (need_hessian) rep(list(matrix(0, p, p)), 4L) else NULL
  hessian_fold_counts <- if (need_hessian) integer(4L) else NULL
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
    score0_i <- as.numeric(-crossprod(X_i, psi_quantile_v2(r_i, tau)) / m_i)
    fold_i <- ((ii - 1L) %% 4L) + 1L

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
      score0 = score0_i,
      fold = fold_i,
      H = H_i,
      H_identity = H_identity_i,
      identity_error = identity_err_i
    )
  }

  cluster_cores <- getOption("jpDMEQR.cluster_cores", 1L)
  out <- if (cluster_cores > 1L && .Platform$OS.type != "windows" && n > 1L) {
    parallel::mclapply(seq_len(n), cluster_worker,
                       mc.cores = min(cluster_cores, n),
                       mc.preschedule = TRUE)
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
    unsmoothed_cluster_scores[ii, ] <- item$score0
    if (need_hessian) {
      H_sum <- H_sum + item$H
      H_identity_sum <- H_identity_sum + item$H_identity
      identity_errors[ii] <- item$identity_error
      f <- item$fold
      hessian_fold_sums[[f]] <- hessian_fold_sums[[f]] + item$H
      hessian_fold_counts[f] <- hessian_fold_counts[f] + 1L
    }
  }

  nuisance_df <- do.call(rbind, nuisance_diagnostics)
  gamma_mat <- do.call(rbind, gamma_list)
  rownames(gamma_mat) <- ci$levels
  colnames(gamma_mat) <- colnames(Z) %||% paste0("z", seq_len(q))

  score <- colMeans(cluster_scores)
  score0 <- colMeans(unsmoothed_cluster_scores)
  hessian <- if (need_hessian) (H_sum / n + t(H_sum / n)) / 2 else NULL
  hessian_identity <- if (need_hessian) (H_identity_sum / n + t(H_identity_sum / n)) / 2 else NULL

  list(
    objective = mean(objective_values),
    score = as.numeric(score),
    score0 = as.numeric(score0),
    hessian = hessian,
    hessian_identity = hessian_identity,
    profile_identity_error = if (need_hessian) max(abs(hessian - hessian_identity)) else NA_real_,
    cluster_identity_errors = identity_errors,
    cluster_scores = cluster_scores,
    unsmoothed_cluster_scores = unsmoothed_cluster_scores,
    hessian_fold_sums = hessian_fold_sums,
    hessian_fold_counts = hessian_fold_counts,
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

# -----------------------------------------------------------------------------
# Round-3 cluster self-normalised first-stage penalty (METHOD_SPECIFICATION_
# ROUND3_AMENDMENT.md sections 1-3).
# -----------------------------------------------------------------------------

# ell_j(b) = [ n^{-1} sum_i { s_ij(b) - sbar_j(b) }^2 ]^{1/2}: the standard
# deviation across INDEPENDENT CLUSTERS of the centred smoothed cluster score
# contributions (n denominator, exactly as specified -- not the n-1 sample sd).
score_loadings_v2 <- function(cluster_scores) {
  S <- as.matrix(cluster_scores)
  if (!nrow(S) || ncol(S) < 1L) stop("cluster_scores must be a non-empty matrix.")
  if (any(!is.finite(S))) {
    stop("Non-finite cluster score contributions in score_loadings_v2.")
  }
  if (nrow(S) == 1L) return(rep(0, ncol(S)))
  S <- S - matrix(colMeans(S), nrow(S), ncol(S), byrow = TRUE)
  sqrt(colMeans(S^2))
}

# Freeze: alpha_lambda,n = 0.10 / log(max(n,3));
#        q_lambda,n      = Phi^{-1}(1 - alpha/(2 p_P));
#        lambda_0,n      = 1.10 q_lambda,n / sqrt(n).
lambda_0_n_v2 <- function(n, p_penalized, safety_constant = 1.10) {
  n <- as.integer(n)
  p_penalized <- as.integer(p_penalized)
  if (n < 2L || p_penalized < 1L) stop("lambda_0_n_v2 requires n>=2 and p_P>=1.")
  alpha <- 0.10 / log(max(n, 3))
  q <- stats::qnorm(1 - alpha / (2 * p_penalized))
  lambda_0 <- safety_constant * q / sqrt(n)
  list(
    alpha = alpha,
    normal_quantile = q,
    safety_constant = safety_constant,
    lambda_0 = lambda_0
  )
}

# Scale-normalised weighted KKT residual R_KKT (round-3 Q1(b)):
#   r_j = |g_j + p_j sign(beta_j)|          j in P, |beta_j| > 1e-10
#       = max(0, |g_j| - p_j)              j in P, |beta_j| <= 1e-10
#       = |g_j|                            j in U (p_j == 0)
#   p_ref = median positive p_j;
#   R_KKT = max( max_{j in P} r_j/p_j, max_{j in U} r_j/p_ref ).
# `penalty` is the EFFECTIVE coordinate penalty vector p_j (zero on U).
kkt_residual_normalized_v2 <- function(beta, score, penalty,
                                       active_tol = 1e-10) {
  beta <- as.numeric(beta)
  score <- as.numeric(score)
  p <- as.numeric(penalty)
  if (length(beta) != length(score) || length(score) != length(p)) {
    stop("beta, score and penalty must have equal length.")
  }
  unpen <- p == 0
  pen <- !unpen
  r <- numeric(length(p))
  active <- pen & abs(beta) > active_tol
  inactive <- pen & !active
  r[active] <- abs(score[active] + p[active] * sign(beta[active]))
  r[inactive] <- pmax(0, abs(score[inactive]) - p[inactive])
  r[unpen] <- abs(score[unpen])
  p_pos <- p[pen & p > 0]
  p_ref <- if (length(p_pos)) stats::median(p_pos) else NA_real_
  if (!is.finite(p_ref) || p_ref <= 0) {
    # No positive penalty anywhere: pure unpenalised KKT is the max |score|.
    return(max(r))
  }
  max(c(
    if (any(pen)) max(r[pen] / p[pen]) else 0,
    if (any(unpen)) max(r[unpen]) / p_ref else 0
  ))
}
