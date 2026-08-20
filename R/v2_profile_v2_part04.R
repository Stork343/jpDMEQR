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

# -----------------------------------------------------------------------------
# Non-truth-based inverse-Hessian CV selection of the Dantzig tolerance
# (theory decision: PILOT_GATE_THEORY_DECISION.md section 3.4).
#
# For each candidate mu on the frozen grid and each of the four deterministic
# cluster folds: solve the Dantzig row on the training-fold Hessian and score
# it on the held-out fold Hessian with the inverse quadratic loss
#     L_f(mu) = 0.5 omega' H_val omega - e_k' omega,
# whose population minimiser is H^{-1} e_k. The selected mu is the LARGEST
# candidate whose mean loss is within one standard error of the minimum.
# No simulation truth, coverage, or bias is used.
select_mu_inverse_hessian_cv_v2 <- function(H,
                                            coordinate,
                                            mu_grid,
                                            fold_sums,
                                            fold_counts,
                                            solver_preference = c("CLARABEL", "ECOS", "SCS"),
                                            solver_opts = list(),
                                            n_folds = 4L) {
  p <- nrow(H)
  H <- (H + t(H)) / 2
  n <- sum(fold_counts)
  e <- numeric(p); e[coordinate] <- 1
  mu_grid <- sort(unique(as.numeric(mu_grid)))
  mu_grid <- mu_grid[mu_grid > 0]

  results <- lapply(mu_grid, function(mu) {
    losses <- numeric(n_folds)
    feasible_folds <- 0L
    for (f in seq_len(n_folds)) {
      n_f <- fold_counts[f]
      if (n_f <= 0L || n_f >= n) next
      H_val_f <- fold_sums[[f]] / n_f
      H_val_f <- (H_val_f + t(H_val_f)) / 2
      H_train_f <- (n * H - n_f * H_val_f) / (n - n_f)
      H_train_f <- (H_train_f + t(H_train_f)) / 2
      sol <- tryCatch(
        solve_dantzig_row_v2(
          H = H_train_f, coordinate = coordinate, mu_grid = c(mu),
          solver_preference = solver_preference, solver_opts = solver_opts
        ),
        error = function(e) list(feasible = FALSE, omega = NULL)
      )
      if (!isTRUE(sol$feasible) || is.null(sol$omega)) next
      feasible_folds <- feasible_folds + 1L
      omega <- as.numeric(sol$omega)
      losses[f] <- 0.5 * drop(crossprod(omega, H_val_f %*% omega)) -
        drop(crossprod(e, omega))
    }
    if (feasible_folds < n_folds) return(NULL)
    list(
      mu = mu,
      loss_mean = mean(losses),
      loss_se = stats::sd(losses) / sqrt(n_folds),
      feasible_folds = feasible_folds
    )
  })
  ok <- results[!vapply(results, is.null, logical(1))]
  if (!length(ok)) {
    return(list(mu = NA_real_, loss_mean = NA_real_, loss_se = NA_real_,
                selected = FALSE, candidates = data.frame()))
  }
  loss_means <- vapply(ok, `[[`, numeric(1), "loss_mean")
  loss_ses <- vapply(ok, `[[`, numeric(1), "loss_se")
  mus <- vapply(ok, `[[`, numeric(1), "mu")
  best <- which.min(loss_means)
  threshold <- loss_means[best] + loss_ses[best]
  in_one_se <- loss_means <= threshold
  # largest feasible multiplier within one SE of the minimum (more regularised)
  selected_idx <- which(in_one_se)[which.max(mus[in_one_se])]
  list(
    mu = mus[selected_idx],
    loss_mean = loss_means[selected_idx],
    loss_se = loss_ses[selected_idx],
    min_loss = loss_means[best],
    selected = TRUE,
    candidates = data.frame(
      mu = mus, loss_mean = loss_means, loss_se = loss_ses,
      in_one_se = in_one_se, stringsAsFactors = FALSE
    )
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

  fold_sums <- fit$components$hessian_fold_sums
  fold_counts <- fit$components$hessian_fold_counts
  use_cv <- !is.null(fold_sums) && length(fold_sums) == 4L &&
    !is.null(fold_counts) && all(fold_counts > 0L)

  for (ii in seq_along(idx)) {
    k <- idx[ii]
    cv <- if (use_cv) {
      select_mu_inverse_hessian_cv_v2(
        H = H, coordinate = k, mu_grid = mu_grid,
        fold_sums = fold_sums, fold_counts = fold_counts,
        solver_preference = solver_preference, solver_opts = solver_opts
      )
    } else {
      list(mu = NA_real_, selected = FALSE,
           candidates = data.frame(mu = NA_real_, loss_mean = NA_real_,
                                   loss_se = NA_real_, in_one_se = NA,
                                   stringsAsFactors = FALSE))
    }

    direction <- if (isTRUE(cv$selected)) {
      solve_dantzig_row_v2(
        H = H, coordinate = k, mu_grid = c(cv$mu),
        solver_preference = solver_preference, solver_opts = solver_opts
      )
    } else {
      # backward-compatible fallback: smallest feasible value on the grid
      solve_dantzig_row_v2(
        H = H, coordinate = k, mu_grid = mu_grid,
        solver_preference = solver_preference, solver_opts = solver_opts
      )
    }
    directions[[ii]] <- direction
    if (!direction$feasible) {
      rows[[ii]] <- data.frame(
        coordinate = names(fit$beta)[k], index = k, feasible = FALSE,
        beta_hat = fit$beta[k], beta_tilde = NA_real_, se = NA_real_,
        se_smoothed = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
        mu = NA_real_, mu_cv_min_loss = cv$min_loss, mu_cv_loss = cv$loss_mean,
        mu_cv_se = cv$loss_se, dantzig_residual = NA_real_,
        omega_l1 = NA_real_, omega_l2 = NA_real_, adjacent_stability = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }

    omega <- direction$omega
    beta_tilde <- as.numeric(fit$beta[k] - sum(omega * fit$components$score))
    # Primary Wald meat: unsmoothed quantile cluster scores (theory decision
    # section 1.1). Bread remains the smoothed effective Hessian via omega.
    score0 <- fit$components$unsmoothed_cluster_scores
    sand0 <- if (!is.null(score0)) {
      cluster_sandwich_coordinate_v2(score0, omega, df_correction = df_correction)
    } else {
      list(se = NA_real_)
    }
    # Diagnostic/sensitivity meat: smoothed cluster scores (section 1.2).
    sandh <- cluster_sandwich_coordinate_v2(
      fit$components$cluster_scores, omega, df_correction = df_correction
    )
    se0 <- sand0$se
    seh <- sandh$se
    se <- if (is.finite(se0)) se0 else seh
    # adjacent-grid direction stability: relative change between consecutive
    # feasible candidate rows (diagnostic only).
    adj <- NA_real_
    cand <- cv$candidates
    if (isTRUE(cv$selected) && nrow(cand) >= 2L) {
      feas_mus <- cand$mu[is.finite(cand$loss_mean)]
      if (length(feas_mus) >= 2L) {
        sols <- lapply(feas_mus, function(mu) {
          s <- tryCatch(
            solve_dantzig_row_v2(H, k, c(mu), solver_preference, solver_opts),
            error = function(e) list(feasible = FALSE)
          )
          if (isTRUE(s$feasible)) as.numeric(s$omega) else NULL
        })
        sols <- sols[!vapply(sols, is.null, logical(1))]
        if (length(sols) >= 2L) {
          adj <- max(vapply(seq_len(length(sols) - 1L), function(j) {
            d <- sols[[j + 1L]] - sols[[j]]
            sqrt(sum(d^2)) / max(sqrt(sum(sols[[j]]^2)), 1e-12)
          }, numeric(1)))
        }
      }
    }
    rows[[ii]] <- data.frame(
      coordinate = names(fit$beta)[k],
      index = k,
      feasible = TRUE,
      beta_hat = as.numeric(fit$beta[k]),
      beta_tilde = beta_tilde,
      se = se,
      se_smoothed = seh,
      ci_lower = beta_tilde - zcrit * se,
      ci_upper = beta_tilde + zcrit * se,
      mu = direction$mu,
      mu_cv_min_loss = cv$min_loss,
      mu_cv_loss = cv$loss_mean,
      mu_cv_se = cv$loss_se,
      dantzig_residual = direction$residual,
      omega_l1 = direction$l1_norm,
      omega_l2 = direction$l2_norm,
      adjacent_stability = adj,
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
