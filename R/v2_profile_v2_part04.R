cluster_sandwich_coordinate_v2 <- function(cluster_scores, omega, df_correction = FALSE) {
  G <- assert_numeric_matrix_v2(cluster_scores, "cluster_scores")
  omega <- as.numeric(omega)
  if (ncol(G) != length(omega)) stop("omega dimension mismatch.")
  n <- nrow(G)
  projected <- as.numeric(G %*% omega)
  projected_centered <- projected - mean(projected)
  sigma2 <- mean(projected_centered^2)
  if (isTRUE(df_correction) && n > 1) sigma2 <- sigma2 * n / (n - 1)
  list(
    sigma2 = sigma2,
    se = sqrt(max(sigma2, 0) / n),
    projected_scores = projected,
    projected_centered = projected_centered,
    n_clusters = n
  )
}

debias_profile_coordinates_v2 <- function(fit,
                                           coordinates,
                                           mu_grid,
                                           ci_level = 0.95,
                                           solver_preference = c("CLARABEL", "ECOS", "SCS"),
                                           solver_opts = list(),
                                           df_correction = FALSE) {
  if (is.null(fit$components$hessian) || is.null(fit$components$score)) {
    stop("fit does not contain profile score/Hessian components.")
  }
  H <- fit$components$hessian
  colnames(H) <- rownames(H) <- names(fit$beta)
  coords <- coordinates
  if (is.character(coords)) {
    idx <- match(coords, names(fit$beta))
    if (anyNA(idx)) stop("Unknown coordinates: ", paste(coords[is.na(idx)], collapse = ", "))
  } else {
    idx <- as.integer(coords)
  }

  alpha <- 1 - ci_level
  zcrit <- stats::qnorm(1 - alpha / 2)
  rows <- vector("list", length(idx))
  directions <- vector("list", length(idx))

  for (ii in seq_along(idx)) {
    k <- idx[ii]
    direction <- solve_dantzig_row_v2(
      H = H, coordinate = k, mu_grid = mu_grid,
      solver_preference = solver_preference,
      solver_opts = solver_opts
    )
    directions[[ii]] <- direction
    if (!direction$feasible) {
      rows[[ii]] <- data.frame(
        coordinate = names(fit$beta)[k], index = k, feasible = FALSE,
        beta_hat = fit$beta[k], beta_tilde = NA_real_, se = NA_real_,
        ci_lower = NA_real_, ci_upper = NA_real_, mu = NA_real_,
        dantzig_residual = NA_real_, omega_l1 = NA_real_, omega_l2 = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }

    omega <- direction$omega
    beta_tilde <- as.numeric(fit$beta[k] - sum(omega * fit$components$score))
    sand <- cluster_sandwich_coordinate_v2(
      fit$components$cluster_scores, omega, df_correction = df_correction
    )
    rows[[ii]] <- data.frame(
      coordinate = names(fit$beta)[k],
      index = k,
      feasible = TRUE,
      beta_hat = as.numeric(fit$beta[k]),
      beta_tilde = beta_tilde,
      se = sand$se,
      ci_lower = beta_tilde - zcrit * sand$se,
      ci_upper = beta_tilde + zcrit * sand$se,
      mu = direction$mu,
      dantzig_residual = direction$residual,
      omega_l1 = direction$l1_norm,
      omega_l2 = direction$l2_norm,
      stringsAsFactors = FALSE
    )
  }

  list(
    table = do.call(rbind, rows),
    directions = directions,
    ci_level = ci_level,
    fit = fit
  )
}

# -----------------------------------------------------------------------------
# Numerical geometry checks
# -----------------------------------------------------------------------------

finite_difference_gradient_v2 <- function(fn, x, eps = 1e-5) {
  x <- as.numeric(x)
  p <- length(x)
  out <- numeric(p)
  for (k in seq_len(p)) {
    step <- numeric(p); step[k] <- eps
    out[k] <- (fn(x + step) - fn(x - step)) / (2 * eps)
  }
  out
}

finite_difference_jacobian_v2 <- function(fn, x, eps = 1e-5) {
  x <- as.numeric(x)
  f0 <- as.numeric(fn(x))
  out <- matrix(NA_real_, nrow = length(f0), ncol = length(x))
  for (k in seq_along(x)) {
    step <- numeric(length(x)); step[k] <- eps
    out[, k] <- (as.numeric(fn(x + step)) - as.numeric(fn(x - step))) / (2 * eps)
  }
  out
}

validate_profile_geometry_v2 <- function(y,
                                         X,
                                         Z,
                                         cluster_id,
                                         beta,
                                         tau,
                                         h,
                                         lambda_gamma = 1,
                                         eps = 1e-5,
                                         nuisance_control = list(
                                           reltol = 1e-13,
                                           maxit = 1000L,
                                           grad_tol = 1e-8
                                         )) {
  analytic <- profile_components_v2(
    y, X, Z, cluster_id, beta, tau, h, lambda_gamma,
    gamma_init = NULL, need_hessian = TRUE,
    nuisance_control = nuisance_control
  )
  obj_fn <- function(b) profile_objective_v2(
    y, X, Z, cluster_id, b, tau, h, lambda_gamma,
    nuisance_control = nuisance_control
  )
  score_fn <- function(b) profile_components_v2(
    y, X, Z, cluster_id, b, tau, h, lambda_gamma,
    need_hessian = FALSE, nuisance_control = nuisance_control
  )$score

  grad_fd <- finite_difference_gradient_v2(obj_fn, beta, eps = eps)
  hess_fd <- finite_difference_jacobian_v2(score_fn, beta, eps = eps)
  hess_fd <- (hess_fd + t(hess_fd)) / 2

  data.frame(
    tau = tau,
    h = h,
    p = ncol(X),
    q = ncol(Z),
    n_clusters = length(unique(cluster_id)),
    gradient_max_error = max(abs(analytic$score - grad_fd)),
    hessian_max_error = max(abs(analytic$hessian - hess_fd)),
    profile_identity_error = analytic$profile_identity_error,
    max_nuisance_gradient = analytic$max_nuisance_gradient,
    nuisance_converged = analytic$nuisance_converged,
    stringsAsFactors = FALSE
  )
}
