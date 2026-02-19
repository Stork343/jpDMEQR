# Epanechnikov kernel density: k(u) = 3/4 (1-u^2) I(|u|<=1)
kernel_epan_pdf <- function(u) {
  out <- numeric(length(u))
  mask <- abs(u) <= 1
  out[mask] <- 0.75 * (1 - u[mask]^2)
  out
}

# CDF of Epanechnikov kernel distribution
kernel_epan_cdf <- function(u) {
  out <- numeric(length(u))
  out[u >= 1] <- 1
  mid <- (u > -1) & (u < 1)
  out[mid] <- 0.5 + 0.75 * (u[mid] - u[mid]^3 / 3)
  out
}

# Primitive of v * k(v): integral v * 3/4*(1-v^2) dv
.kernel_epan_m1_primitive <- function(v) {
  (3 / 8) * v^2 - (3 / 16) * v^4
}

# Smoothed check loss rho_h(u) = E[rho_tau(u - hV)], V~Epanechnikov
rho_smooth <- function(u, tau, h) {
  if (h <= 0) stop("h must be positive.")
  x <- u / h

  # Region split by v <= x and v > x on support [-1, 1]
  upper1 <- pmin(1, pmax(-1, x))
  lower2 <- upper1

  A1 <- kernel_epan_cdf(upper1) - kernel_epan_cdf(rep(-1, length(u)))
  A2 <- kernel_epan_cdf(rep(1, length(u))) - kernel_epan_cdf(lower2)

  M_minus1 <- .kernel_epan_m1_primitive(-1)
  M_plus1 <- .kernel_epan_m1_primitive(1)
  B1 <- .kernel_epan_m1_primitive(upper1) - M_minus1
  B2 <- M_plus1 - .kernel_epan_m1_primitive(lower2)

  tau * (u * A1 - h * B1) + (tau - 1) * (u * A2 - h * B2)
}

# First derivative psi_h(u)
psi_smooth <- function(u, tau, h) {
  if (h <= 0) stop("h must be positive.")
  tau - kernel_epan_cdf(-u / h)
}

# Second derivative phi_h(u)
phi_smooth <- function(u, tau, h) {
  if (h <= 0) stop("h must be positive.")
  kernel_epan_pdf(u / h) / h
}

# ==============================================================================
# 2) Cluster-adjusted IQR-SIS (multi-quantile)
# ==============================================================================

#' Cluster-Adjusted IQR-SIS Screening
#'
#' Iterative screening for clustered quantile-regression data using within-cluster
#' median centering and marginal multi-quantile scores.
#'
#' @param y Numeric response vector of length `N`.
#' @param X Numeric design matrix of dimension `N x p`.
#' @param cluster_id Cluster labels of length `N`.
#' @param tau_grid Quantile levels used to build marginal screening scores.
#' @param tau_target Target quantile used in the residual update step.
#' @param d_target Target screening size.
#' @param T_iter Number of iterative screening rounds.
#' @param cumulative_targets Optional cumulative screening sizes per round.
#' @param verbose Logical; print progress if `TRUE`.
#'
#' @return Integer vector of selected variable indices.
#' @export
CA_IQR_SIS <- function(y,
                       X,
                       cluster_id,
                       tau_grid = c(0.25, 0.5, 0.75),
                       tau_target = 0.5,
                       d_target,
                       T_iter = 2,
                       cumulative_targets = NULL,
                       verbose = FALSE) {
  y <- as.numeric(y)
  X <- as.matrix(X)
  cluster_id <- as.vector(cluster_id)

  N <- length(y)
  p <- ncol(X)
  if (length(cluster_id) != N) stop("cluster_id length mismatch.")
  if (missing(d_target) || is.null(d_target)) stop("d_target must be provided.")

  tau_grid <- sort(unique(as.numeric(tau_grid)))
  if (any(tau_grid <= 0 | tau_grid >= 1)) stop("tau_grid must be in (0,1).")
  if (tau_target <= 0 || tau_target >= 1) stop("tau_target must be in (0,1).")

  d_target <- max(1L, min(as.integer(d_target), p))
  T_iter <- max(1L, as.integer(T_iter))

  if (is.null(cumulative_targets)) {
    cumulative_targets <- round(seq(d_target / T_iter, d_target, length.out = T_iter))
  }
  cumulative_targets <- as.integer(cumulative_targets)
  if (length(cumulative_targets) != T_iter) {
    stop("cumulative_targets must have length T_iter.")
  }
  cumulative_targets <- pmax(0L, pmin(d_target, cumulative_targets))
  cumulative_targets <- cummax(cumulative_targets)

  active_set <- integer(0)
  residuals <- y

  for (t in seq_len(T_iter)) {
    # 1) within-cluster median centering
    cluster_medians <- tapply(residuals, cluster_id, median)
    mu_hat <- as.numeric(cluster_medians[as.character(cluster_id)])
    r_tilde <- residuals - mu_hat

    # 2) marginal multi-quantile score I_k = max_l |beta_k(tau_l)|
    candidates <- setdiff(seq_len(p), active_set)
    if (length(candidates) == 0) break

    scores <- rep(0, p)
    for (j in candidates) {
      xj <- X[, j]
      tau_scores <- numeric(length(tau_grid))

      for (k in seq_along(tau_grid)) {
        tau_k <- tau_grid[k]
        fit <- tryCatch(
          quantreg::rq(r_tilde ~ xj, tau = tau_k, method = "fn"),
          error = function(e) NULL
        )
        if (!is.null(fit)) {
          coef_k <- coef(fit)
          if (length(coef_k) >= 2 && is.finite(coef_k[2])) {
            tau_scores[k] <- abs(coef_k[2])
          }
        }
      }
      scores[j] <- max(tau_scores)
    }

    # 3) add variables to reach cumulative target of round t
    n_to_add <- max(0L, cumulative_targets[t] - length(active_set))
    if (n_to_add > 0L) {
      ranked <- candidates[order(scores[candidates], decreasing = TRUE)]
      new_active <- ranked[seq_len(min(n_to_add, length(ranked)))]
      active_set <- sort(unique(c(active_set, new_active)))
    }

    # 4) residual update at target quantile
    if (length(active_set) > 0) {
      X_sub <- X[, active_set, drop = FALSE]
      fit_update <- tryCatch(
        quantreg::rq(y ~ ., data = data.frame(y = y, X_sub), tau = tau_target, method = "fn"),
        error = function(e) NULL
      )
      if (!is.null(fit_update)) {
        residuals <- as.numeric(fit_update$residuals)
      }
    }

    if (verbose) {
      cat(sprintf("[CA-IQR-SIS] round %d/%d: active=%d\n", t, T_iter, length(active_set)))
    }
  }

  active_set
}

# ==============================================================================
# 3) Joint penalised estimation via smoothed BCD
# ==============================================================================

soft_threshold <- function(x, lam) {
  sign(x) * pmax(0, abs(x) - lam)
}

# Theory-scale CLIME tuning: sqrt(log p / (N h)) + h^2
#' Theory-Scale CLIME Tuning Level
#'
#' Computes the default CLIME tuning level
#' `sqrt(log(p)/(N h)) + h^2`.
#'
#' @param p Number of fixed-effect covariates.
#' @param N Total sample size.
#' @param h Smoothing bandwidth.
#'
#' @return Numeric scalar tuning level.
#' @export
default_lambda_clime <- function(p, N, h) {
  if (h <= 0) stop("h must be positive.")
  sqrt(log(max(2, p)) / (N * h)) + h^2
}

# Gamma subproblem: warm-start BFGS (line-search handled by optim)
solve_gamma_subproblem <- function(r_fixed,
                                   Z_i,
                                   tau,
                                   h,
                                   lambda2,
                                   N,
                                   gamma_init = NULL,
                                   maxit = 100,
                                   reltol = 1e-8) {
  q <- ncol(Z_i)
  if (is.null(gamma_init)) gamma_init <- rep(0, q)

  obj_fn <- function(gamma_i) {
    res <- r_fixed - as.numeric(Z_i %*% gamma_i)
    loss <- sum(rho_smooth(res, tau, h)) / N
    pen <- (lambda2 / (2 * N)) * sum(gamma_i^2)
    loss + pen
  }

  grad_fn <- function(gamma_i) {
    res <- r_fixed - as.numeric(Z_i %*% gamma_i)
    psi <- psi_smooth(res, tau, h)
    as.numeric(-(t(Z_i) %*% psi) / N + (lambda2 / N) * gamma_i)
  }

  opt <- optim(
    par = gamma_init,
    fn = obj_fn,
    gr = grad_fn,
    method = "BFGS",
    control = list(maxit = maxit, reltol = reltol)
  )

  list(
    gamma = as.numeric(opt$par),
    value = opt$value,
    convergence = opt$convergence,
    fn_counts = as.numeric(opt$counts["function"]),
    gr_counts = as.numeric(opt$counts["gradient"])
  )
}

# Beta subproblem: coordinate descent with KKT check
solve_beta_subproblem <- function(y_adj,
                                  X,
                                  beta_init,
                                  tau,
                                  h,
                                  lambda1,
                                  N,
                                  max_iter = 200,
                                  tol = 1e-6,
                                  kkt_tol = 1e-5) {
  X <- as.matrix(X)
  p <- ncol(X)
  beta <- as.numeric(beta_init)
  if (length(beta) != p) stop("beta_init length mismatch.")

  # ISTA with a global Lipschitz constant for stable updates.
  max_phi <- 0.75 / h
  xtx_n <- crossprod(X) / N
  lmax <- .largest_eigenvalue_power(xtx_n)
  L_global <- max(max_phi * lmax, 1e-6)

  r <- as.numeric(y_adj - X %*% beta)
  converged <- FALSE

  for (iter in seq_len(max_iter)) {
    beta_old <- beta
    psi <- psi_smooth(r, tau, h)
    grad <- as.numeric(-crossprod(X, psi) / N)
    beta <- soft_threshold(beta - grad / L_global, lambda1 / L_global)
    r <- as.numeric(y_adj - X %*% beta)

    if (max(abs(beta - beta_old)) < tol) {
      psi_now <- psi_smooth(r, tau, h)
      grad <- as.numeric(-crossprod(X, psi_now) / N)

      active <- abs(beta) > 0
      kkt_violation <- numeric(p)
      if (any(active)) {
        kkt_violation[active] <- abs(grad[active] + lambda1 * sign(beta[active]))
      }
      if (any(!active)) {
        kkt_violation[!active] <- pmax(0, abs(grad[!active]) - lambda1)
      }

      if (max(kkt_violation) < kkt_tol) {
        converged <- TRUE
        break
      }
    }
  }

  list(beta = beta, residual = r, converged = converged)
}

#' Joint Penalised Mixed-Effects Quantile Regression
#'
#' Fits the smoothed objective using block coordinate descent with
#' an `L1` penalty on fixed effects and ridge penalty on random effects.
#'
#' @param y Numeric response vector.
#' @param X Fixed-effect design matrix.
#' @param Z Random-effect design matrix.
#' @param cluster_id Cluster labels.
#' @param tau Target quantile in `(0, 1)`.
#' @param h Smoothing bandwidth.
#' @param lambda1 `L1` penalty level for fixed effects.
#' @param lambda2 Ridge penalty level for random effects.
#' @param max_iter Maximum outer BCD iterations.
#' @param tol Outer-loop convergence tolerance for parameter changes.
#' @param obj_tol Relative objective-change tolerance.
#' @param initial_beta Optional initialization for `beta`.
#' @param beta_cd_max_iter Maximum inner iterations for beta subproblem.
#' @param beta_cd_tol Convergence tolerance for beta subproblem.
#' @param verbose Logical; print progress if `TRUE`.
#'
#' @return A list containing penalised estimates, diagnostics, and convergence info.
#' @export
JP_DME_QR <- function(y,
                      X,
                      Z,
                      cluster_id,
                      tau,
                      h,
                      lambda1,
                      lambda2,
                      max_iter = 100,
                      tol = 1e-4,
                      obj_tol = 1e-6,
                      initial_beta = NULL,
                      beta_cd_max_iter = 200,
                      beta_cd_tol = 1e-6,
                      verbose = FALSE) {
  y <- as.numeric(y)
  X <- as.matrix(X)
  Z <- as.matrix(Z)

  N <- length(y)
  p <- ncol(X)
  q <- ncol(Z)
  if (nrow(X) != N || nrow(Z) != N || length(cluster_id) != N) {
    stop("Dimension mismatch among y, X, Z, cluster_id.")
  }

  cluster_fac <- factor(cluster_id)
  cluster_levels <- levels(cluster_fac)
  n_clusters <- nlevels(cluster_fac)
  cluster_rows <- split(seq_len(N), cluster_fac)

  if (is.null(initial_beta)) initial_beta <- rep(0, p)
  beta <- as.numeric(initial_beta)
  gamma_mat <- matrix(0, nrow = n_clusters, ncol = q)
  rownames(gamma_mat) <- cluster_levels

  joint_objective <- function(beta_now, gamma_now) {
    Z_gamma_now <- numeric(N)
    for (i in seq_len(n_clusters)) {
      idx <- cluster_rows[[i]]
      Z_gamma_now[idx] <- as.numeric(Z[idx, , drop = FALSE] %*% gamma_now[i, ])
    }
    res_now <- y - as.numeric(X %*% beta_now) - Z_gamma_now
    mean_loss <- sum(rho_smooth(res_now, tau, h)) / N
    l1_pen <- lambda1 * sum(abs(beta_now))
    ridge_pen <- (lambda2 / (2 * N)) * sum(gamma_now^2)
    mean_loss + l1_pen + ridge_pen
  }

  obj_prev <- joint_objective(beta, gamma_mat)
  objective_history <- rep(NA_real_, max_iter + 1L)
  objective_history[1L] <- obj_prev
  avg_bfgs_steps_history <- rep(NA_real_, max_iter)

  converged <- FALSE
  n_iter <- max_iter

  for (k in seq_len(max_iter)) {
    beta_old <- beta

    # Step 1: beta update given gamma
    Z_gamma <- numeric(N)
    for (i in seq_len(n_clusters)) {
      idx <- cluster_rows[[i]]
      Z_gamma[idx] <- as.numeric(Z[idx, , drop = FALSE] %*% gamma_mat[i, ])
    }
    y_adj <- y - Z_gamma

    beta_fit <- solve_beta_subproblem(
      y_adj = y_adj,
      X = X,
      beta_init = beta,
      tau = tau,
      h = h,
      lambda1 = lambda1,
      N = N,
      max_iter = beta_cd_max_iter,
      tol = beta_cd_tol
    )
    beta <- beta_fit$beta

    # Step 2: gamma update given beta (warm start)
    r_fixed <- y - as.numeric(X %*% beta)
    max_gamma_diff <- 0
    bfgs_fn_counts <- numeric(n_clusters)

    for (i in seq_len(n_clusters)) {
      idx <- cluster_rows[[i]]
      Zi <- Z[idx, , drop = FALSE]
      ri <- r_fixed[idx]

      gamma_fit <- solve_gamma_subproblem(
        r_fixed = ri,
        Z_i = Zi,
        tau = tau,
        h = h,
        lambda2 = lambda2,
        N = N,
        gamma_init = gamma_mat[i, ]
      )
      gamma_new <- gamma_fit$gamma

      diff_i <- sqrt(sum((gamma_new - gamma_mat[i, ])^2))
      max_gamma_diff <- max(max_gamma_diff, diff_i)
      gamma_mat[i, ] <- gamma_new
      bfgs_fn_counts[i] <- gamma_fit$fn_counts
    }
    avg_bfgs_steps_history[k] <- mean(bfgs_fn_counts, na.rm = TRUE)

    beta_diff <- max(abs(beta - beta_old))
    obj_now <- joint_objective(beta, gamma_mat)
    objective_history[k + 1L] <- obj_now
    obj_rel_change <- abs(obj_now - obj_prev) / (1 + abs(obj_prev))

    if (verbose) {
      cat(sprintf(
        "[BCD] iter %d: beta_inf_diff=%.3e, gamma_l2_diff=%.3e, rel_obj=%.3e\n",
        k, beta_diff, max_gamma_diff, obj_rel_change
      ))
    }

    if (beta_diff < tol && max_gamma_diff < tol && obj_rel_change < obj_tol) {
      converged <- TRUE
      n_iter <- k
      break
    }
    obj_prev <- obj_now
  }

  if (!converged) {
    n_iter <- max_iter
  }

  objective_history <- objective_history[seq_len(n_iter + 1L)]
  avg_bfgs_steps_history <- avg_bfgs_steps_history[seq_len(n_iter)]

  list(
    beta = beta,
    gamma = gamma_mat,
    cluster_levels = cluster_levels,
    converged = converged,
    iterations = n_iter,
    objective_history = objective_history,
    average_bfgs_steps = if (length(avg_bfgs_steps_history) == 0) NA_real_ else mean(avg_bfgs_steps_history, na.rm = TRUE),
    average_bfgs_steps_history = avg_bfgs_steps_history
  )
}

# ==============================================================================
# 4) Debiased inference (Schur complement + CLIME)
# ==============================================================================

.align_gamma_matrix <- function(gamma_hat, cluster_id) {
  cluster_fac <- factor(cluster_id)
  lev <- levels(cluster_fac)

  if (is.vector(gamma_hat)) {
    gamma_hat <- matrix(gamma_hat, ncol = 1)
  }
  gamma_hat <- as.matrix(gamma_hat)

  if (!is.null(rownames(gamma_hat))) {
    if (!all(lev %in% rownames(gamma_hat))) {
      stop("gamma_hat rownames do not cover all cluster levels.")
    }
    return(gamma_hat[lev, , drop = FALSE])
  }

  if (nrow(gamma_hat) != length(lev)) {
    stop("gamma_hat rows must match number of clusters.")
  }
  rownames(gamma_hat) <- lev
  gamma_hat
}

Compute_Hessian_Score <- function(y,
                                  X,
                                  Z,
                                  cluster_id,
                                  beta_hat,
                                  gamma_hat,
                                  tau,
                                  h,
                                  lambda2) {
  y <- as.numeric(y)
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  N <- length(y)
  p <- ncol(X)
  q <- ncol(Z)

  cluster_fac <- factor(cluster_id)
  cluster_levels <- levels(cluster_fac)
  n_clusters <- nlevels(cluster_fac)
  cluster_rows <- split(seq_len(N), cluster_fac)

  gamma_mat <- .align_gamma_matrix(gamma_hat, cluster_id)

  Z_gamma <- numeric(N)
  for (i in seq_len(n_clusters)) {
    idx <- cluster_rows[[i]]
    Z_gamma[idx] <- as.numeric(Z[idx, , drop = FALSE] %*% gamma_mat[i, ])
  }

  res <- y - as.numeric(X %*% beta_hat) - Z_gamma
  psi_vals <- psi_smooth(res, tau, h)
  phi_vals <- phi_smooth(res, tau, h)

  H_bb <- (t(X) %*% (X * phi_vals)) / N
  H_bg_H_gg_inv_H_gb <- matrix(0, p, p)
  X_tilde <- matrix(0, N, p)
  H_gg_condition <- rep(NA_real_, n_clusters)

  for (i in seq_len(n_clusters)) {
    idx <- cluster_rows[[i]]
    Xi <- X[idx, , drop = FALSE]
    Zi <- Z[idx, , drop = FALSE]
    wi <- phi_vals[idx]

    H_gg_i <- (t(Zi) %*% (Zi * wi)) / N + diag(lambda2 / N, q)
    H_bg_i <- (t(Xi) %*% (Zi * wi)) / N
    H_gg_condition[i] <- tryCatch(kappa(H_gg_i), error = function(e) Inf)

    inv_H_gg_i <- tryCatch(
      solve(H_gg_i),
      error = function(e) solve(H_gg_i + diag(1e-8, q))
    )

    H_gb_i <- t(H_bg_i)
    H_bg_H_gg_inv_H_gb <- H_bg_H_gg_inv_H_gb + H_bg_i %*% inv_H_gg_i %*% H_gb_i

    correction <- Zi %*% inv_H_gg_i %*% H_gb_i
    X_tilde[idx, ] <- Xi - correction
  }

  H_eff <- H_bb - H_bg_H_gg_inv_H_gb
  H_eff_condition <- tryCatch(kappa(H_eff), error = function(e) Inf)
  g_hat <- as.numeric(-crossprod(X_tilde, psi_vals) / N)

  U_i <- matrix(0, nrow = n_clusters, ncol = p)
  for (i in seq_len(n_clusters)) {
    idx <- cluster_rows[[i]]
    U_i[i, ] <- as.numeric(crossprod(X_tilde[idx, , drop = FALSE], psi_vals[idx]))
  }

  J_hat <- crossprod(U_i)
  # Cluster-robust "meat" for g_hat = N^{-1} sum_i U_i.
  Sigma_xi_hat <- J_hat / (N^2)

  list(
    H_eff = H_eff,
    g_hat = g_hat,
    J_hat = J_hat,
    Sigma_xi_hat = Sigma_xi_hat,
    X_tilde = X_tilde,
    residuals = res,
    cluster_levels = cluster_levels,
    H_gg_condition = H_gg_condition,
    H_eff_condition = H_eff_condition
  )
}

# Largest eigenvalue approximation for Lipschitz constant in ISTA
.largest_eigenvalue_power <- function(A, max_iter = 200, tol = 1e-7) {
  p <- ncol(A)
  v <- rep(1 / sqrt(p), p)
  lambda_old <- 0
  for (k in seq_len(max_iter)) {
    Av <- as.numeric(A %*% v)
    nrm <- sqrt(sum(Av^2))
    if (nrm < 1e-12) return(0)
    v <- Av / nrm
    lambda_new <- as.numeric(crossprod(v, A %*% v))
    if (abs(lambda_new - lambda_old) < tol * max(1, abs(lambda_new))) {
      return(max(lambda_new, 0))
    }
    lambda_old <- lambda_new
  }
  max(lambda_old, 0)
}

# Solve min_w 0.5||H w - b||_2^2 + alpha ||w||_1 by ISTA/FISTA
.solve_lasso_ista <- function(H,
                              b,
                              alpha,
                              w_init = NULL,
                              max_iter = 500,
                              tol = 1e-6) {
  p <- ncol(H)
  if (is.null(w_init)) w_init <- rep(0, p)

  HtH <- crossprod(H)
  L <- .largest_eigenvalue_power(HtH)
  L <- max(L, 1e-8)

  w <- as.numeric(w_init)
  z <- w
  t_now <- 1
  Htb <- as.numeric(crossprod(H, b))

  for (k in seq_len(max_iter)) {
    w_old <- w
    grad <- as.numeric(HtH %*% z - Htb)
    w <- soft_threshold(z - grad / L, alpha / L)

    t_next <- (1 + sqrt(1 + 4 * t_now^2)) / 2
    z <- w + ((t_now - 1) / t_next) * (w - w_old)
    t_now <- t_next

    if (max(abs(w - w_old)) < tol) break
  }

  w
}

# ADMM fallback for exact CLIME column:
# min ||w||_1 s.t. ||H w - e_j||_inf <= lambda
solve_clime_column_admm <- function(H_eff,
                                    target_col,
                                    lambda_clime,
                                    rho = 1,
                                    max_iter = 2000,
                                    tol = 1e-4,
                                    ista_max_iter = 500,
                                    ista_tol = 1e-6) {
  p <- ncol(H_eff)
  e_j <- numeric(p)
  e_j[target_col] <- 1

  w <- rep(0, p)
  r <- rep(0, p)
  u <- rep(0, p)

  for (k in seq_len(max_iter)) {
    # w-update: lasso subproblem
    b <- e_j + r - u
    w <- .solve_lasso_ista(
      H = H_eff,
      b = b,
      alpha = 1 / rho,
      w_init = w,
      max_iter = ista_max_iter,
      tol = ista_tol
    )

    Hw_minus_e <- as.numeric(H_eff %*% w - e_j)
    r_old <- r

    # r-update: projection onto infinity-norm box
    r <- pmin(lambda_clime, pmax(-lambda_clime, Hw_minus_e + u))

    # dual update
    u <- u + Hw_minus_e - r

    primal_res <- max(abs(Hw_minus_e - r))
    dual_res <- max(abs(r - r_old))
    if (primal_res < tol && dual_res < tol) break
  }

  w
}

# Exact CLIME column solver: CVXR if available, otherwise ADMM fallback
solve_clime_column <- function(H_eff,
                               target_col,
                               lambda_clime,
                               solver_order = c("ECOS", "SCS"),
                               method = c("auto", "cvxr", "admm"),
                               admm_control = list(),
                               verbose = FALSE) {
  method <- match.arg(method)
  if (method == "auto") {
    method <- if (requireNamespace("CVXR", quietly = TRUE)) "cvxr" else "admm"
  }

  if (method == "cvxr") {
    if (!requireNamespace("CVXR", quietly = TRUE)) {
      stop("method='cvxr' requested but package 'CVXR' is not installed.")
    }
    cvxr_ns <- asNamespace("CVXR")
    cvxr_variable <- get("Variable", envir = cvxr_ns)
    cvxr_minimize <- get("Minimize", envir = cvxr_ns)
    cvxr_p_norm <- get("p_norm", envir = cvxr_ns)
    cvxr_problem <- get("Problem", envir = cvxr_ns)
    cvxr_solve <- get("solve", envir = cvxr_ns)

    p <- ncol(H_eff)
    e_j <- numeric(p)
    e_j[target_col] <- 1

    w <- cvxr_variable(p)
    obj <- cvxr_minimize(cvxr_p_norm(w, 1))
    cons <- list(cvxr_p_norm(H_eff %*% w - e_j, "inf") <= lambda_clime)
    prob <- cvxr_problem(obj, cons)

    for (solver_name in solver_order) {
      fit <- tryCatch(
        cvxr_solve(prob, solver = solver_name, verbose = verbose),
        error = function(e) NULL
      )

      if (!is.null(fit) && fit$status %in% c("optimal", "optimal_inaccurate")) {
        return(as.numeric(fit$getValue(w)))
      }
    }
    stop(sprintf("CLIME column %d optimization failed for all CVXR solvers.", target_col))
  }

  # ADMM fallback
  ctrl <- modifyList(
    list(rho = 1, max_iter = 2000, tol = 1e-4, ista_max_iter = 500, ista_tol = 1e-6),
    admm_control
  )
  solve_clime_column_admm(
    H_eff = H_eff,
    target_col = target_col,
    lambda_clime = lambda_clime,
    rho = ctrl$rho,
    max_iter = ctrl$max_iter,
    tol = ctrl$tol,
    ista_max_iter = ctrl$ista_max_iter,
    ista_tol = ctrl$ista_tol
  )
}

estimate_precision_matrix_clime <- function(H_eff,
                                            lambda_clime,
                                            solver_order = c("ECOS", "SCS"),
                                            clime_method = c("auto", "cvxr", "admm"),
                                            admm_control = list(),
                                            verbose = FALSE) {
  clime_method <- match.arg(clime_method)
  p <- ncol(H_eff)
  Omega <- matrix(0, p, p)

  for (j in seq_len(p)) {
    if (verbose) cat(sprintf("[CLIME] Solving column %d/%d\n", j, p))
    Omega[, j] <- solve_clime_column(
      H_eff = H_eff,
      target_col = j,
      lambda_clime = lambda_clime,
      solver_order = solver_order,
      method = clime_method,
      admm_control = admm_control,
      verbose = FALSE
    )
  }

  # Symmetrization
  0.5 * (Omega + t(Omega))
}

Debias_Inference.CLIME <- function(y,
                                   X,
                                   Z,
                                   cluster_id,
                                   beta_hat,
                                   gamma_hat,
                                   tau,
                                   h,
                                   lambda2,
                                   lambda_clime = NULL,
                                   clime_solver_order = c("ECOS", "SCS"),
                                   clime_method = c("auto", "cvxr", "admm"),
                                   clime_admm_control = list(),
                                   verbose = FALSE) {
  N <- length(y)
  p <- ncol(as.matrix(X))
  clime_method <- match.arg(clime_method)

  if (is.null(lambda_clime)) {
    lambda_clime <- default_lambda_clime(p = p, N = N, h = h)
  }

  comp <- Compute_Hessian_Score(
    y = y,
    X = X,
    Z = Z,
    cluster_id = cluster_id,
    beta_hat = beta_hat,
    gamma_hat = gamma_hat,
    tau = tau,
    h = h,
    lambda2 = lambda2
  )

  Omega_hat <- estimate_precision_matrix_clime(
    H_eff = comp$H_eff,
    lambda_clime = lambda_clime,
    solver_order = clime_solver_order,
    clime_method = clime_method,
    admm_control = clime_admm_control,
    verbose = verbose
  )

  beta_tilde <- as.numeric(beta_hat - Omega_hat %*% comp$g_hat)

  # Cluster-robust sandwich variance with N^{-1} score scaling.
  V <- Omega_hat %*% comp$Sigma_xi_hat %*% t(Omega_hat)
  se <- sqrt(pmax(diag(V), 0))

  list(
    beta_tilde = beta_tilde,
    se = se,
    Omega = Omega_hat,
    H_eff = comp$H_eff,
    g_hat = comp$g_hat,
    Sigma_xi = comp$Sigma_xi_hat,
    lambda_clime = lambda_clime,
    variance_matrix = V,
    H_eff_condition = comp$H_eff_condition,
    H_gg_condition = comp$H_gg_condition
  )
}

Debias_Inference.Ridge <- function(y,
                                   X,
                                   Z,
                                   cluster_id,
                                   beta_hat,
                                   gamma_hat,
                                   tau,
                                   h,
                                   lambda2,
                                   ridge_const = 1e-4) {
  N <- length(y)
  p <- ncol(as.matrix(X))

  comp <- Compute_Hessian_Score(
    y = y,
    X = X,
    Z = Z,
    cluster_id = cluster_id,
    beta_hat = beta_hat,
    gamma_hat = gamma_hat,
    tau = tau,
    h = h,
    lambda2 = lambda2
  )

  H_reg <- comp$H_eff + diag(ridge_const, p)
  Omega_hat <- tryCatch(
    solve(H_reg),
    error = function(e) solve(H_reg + diag(1e-3, p))
  )

  beta_tilde <- as.numeric(beta_hat - Omega_hat %*% comp$g_hat)

  V <- Omega_hat %*% comp$Sigma_xi_hat %*% t(Omega_hat)
  se <- sqrt(pmax(diag(V), 0))

  list(
    beta_tilde = beta_tilde,
    se = se,
    Omega = Omega_hat,
    H_eff = comp$H_eff,
    g_hat = comp$g_hat,
    Sigma_xi = comp$Sigma_xi_hat,
    variance_matrix = V,
    H_eff_condition = comp$H_eff_condition,
    H_gg_condition = comp$H_gg_condition
  )
}

# Unified API used by simulations and real-data script
#' Debiased Inference for Penalised Mixed-Effects QR
#'
#' Applies Schur-complement adjustment and precision-matrix estimation
#' (CLIME or ridge inverse) to obtain debiased estimates and standard errors.
#'
#' @param y Numeric response vector.
#' @param X Fixed-effect design matrix.
#' @param Z Random-effect design matrix.
#' @param cluster_id Cluster labels.
#' @param beta_hat Penalised fixed-effect estimate.
#' @param gamma_hat Penalised random-effect estimate matrix by cluster.
#' @param tau Target quantile.
#' @param h Smoothing bandwidth.
#' @param lambda2 Ridge penalty level used in fitting.
#' @param lambda_clime Optional CLIME tuning level.
#' @param method Debiasing method: `"CLIME"` or `"RIDGE"`.
#' @param ... Additional arguments passed to the chosen method.
#'
#' @return A list with debiased estimates, standard errors, and variance diagnostics.
#' @export
Debias_Inference <- function(y,
                             X,
                             Z,
                             cluster_id,
                             beta_hat,
                             gamma_hat,
                             tau,
                             h,
                             lambda2,
                             lambda_clime = NULL,
                             method = c("CLIME", "RIDGE"),
                             ...) {
  method <- match.arg(method)

  if (method == "CLIME") {
    Debias_Inference.CLIME(
      y = y,
      X = X,
      Z = Z,
      cluster_id = cluster_id,
      beta_hat = beta_hat,
      gamma_hat = gamma_hat,
      tau = tau,
      h = h,
      lambda2 = lambda2,
      lambda_clime = lambda_clime,
      ...
    )
  } else {
    Debias_Inference.Ridge(
      y = y,
      X = X,
      Z = Z,
      cluster_id = cluster_id,
      beta_hat = beta_hat,
      gamma_hat = gamma_hat,
      tau = tau,
      h = h,
      lambda2 = lambda2,
      ...
    )
  }
}
