# -----------------------------------------------------------------------------
# Oracle variance ladder (ROUND-2 adjudication, item A1).
#
# For the debiased estimator at coordinate k define
#     T_k = sqrt(n) ( beta_tilde_k - beta*_k ),
#     L_k = -n^{-1/2} sum_i omega_k^{pop}' g_i^*,      (oracle influence at target)
#     R_k = T_k - L_k.
# We report Var(T), Var(L), Var(R) and 2Cov(L,R) across replications, and the
# four-level score ladder that localises where the finite-sample SE shortfall
# arises:
#   level 1  population first-order meat      = omega^{pop}' Sigma_pop omega^{pop}
#   level 2  target-score meat (pop direction) = omega^{pop}' Sigma(gastar) omega^{pop}
#   level 3  fitted-score meat (pop direction) = omega^{pop}' Sigma(gfitted0) omega^{pop}
#   level 4  fitted-score meat (sample exact)  = omega^{exact}' Sigma(gfitted0) omega^{exact}
# where Sigma_pop is the population unsmoothed score covariance from the POP-H
# asset, g_i^* uses the UNsmoothed score at the profile-optimal residual for
# the target beta*, g_i^(0) uses the fitted residuals, and omega^{exact} is the
# sample exact inverse row H^{-1} e_k. These are diagnostics only; none is used
# to select mu or to tune the primary interval.
# -----------------------------------------------------------------------------

# Target residual scores g_i^* = -m_i^{-1} X_i' psi_tau( r_i(beta_star) ).
# r_i(beta*) = y_i - X_i beta* - Z_i gamma_hat(beta*), with gamma minimising the
# (smoothed) cluster criterion at beta*.
target_residual_scores_v2 <- function(dat, beta_target, tau, lambda_gamma,
                                      h = 0, nuisance_control = list(reltol = 1e-10,
                                                                      maxit = 300L,
                                                                      grad_tol = 1e-7)) {
  X <- dat$X; Z <- dat$Z; y <- dat$y; cid <- dat$cluster_id
  ci <- cluster_index_v2(cid)
  Lambda <- make_lambda_matrix_v2(lambda_gamma, ncol(Z))
  G <- matrix(NA_real_, nrow = ci$n, ncol = ncol(X),
              dimnames = list(ci$levels, colnames(X)))
  for (ii in seq_len(ci$n)) {
    idx <- ci$rows[[ii]]
    X_i <- X[idx, , drop = FALSE]
    Z_i <- Z[idx, , drop = FALSE]
    m_i <- length(idx)
    gam <- profile_cluster_gamma_v2(
      y = y[idx], X_i = X_i, Z_i = Z_i, beta = as.numeric(beta_target),
      tau = tau, h = h, Lambda = Lambda, gamma_init = NULL,
      control = modifyList(list(reltol = 1e-8, maxit = 100L, grad_tol = 1e-6),
                           nuisance_control)
    )$gamma
    r_i <- y[idx] - as.numeric(X_i %*% beta_target) - as.numeric(Z_i %*% gam)
    G[ii, ] <- as.numeric(-crossprod(X_i, psi_quantile_v2(r_i, tau)) / m_i)
  }
  G
}

# One coordinate's ladder. Returns meat-per-cluster levels and L_k.
# Gs (target unsmoothed scores at beta*) is precomputed ONCE per replication by
# the caller and shared across coordinates to avoid repeating the full profile
# pass.
variance_ladder_one_v2 <- function(dat, fit, beta_target, population_asset,
                                   coordinate, tau, lambda_gamma,
                                   Gs = NULL,
                                   h = fit$fit_object$components$h %||% 0,
                                   nuisance_control = NULL) {
  fit_names <- names(fit$beta_hat)
  k <- match(coordinate, fit_names)
  e <- numeric(length(fit_names)); e[k] <- 1
  Hf <- fit$fit_object$components$hessian
  G0 <- fit$fit_object$components$unsmoothed_cluster_scores
  n <- nrow(G0)
  out <- list(L_k = NA_real_, meat1 = NA_real_, meat2 = NA_real_,
              meat3 = NA_real_, meat4 = NA_real_)
  # population direction + population meat (asset)
  if (is.null(population_asset) || is.null(population_asset$directions) ||
      !coordinate %in% names(population_asset$directions)) return(out)
  pop <- population_asset$directions[[coordinate]]$omega
  pop_full <- setNames(rep(0, length(fit_names)), fit_names)
  common <- intersect(names(pop), fit_names)
  pop_full[common] <- pop[common]
  Sigma_pop <- population_asset$Sigma0_population
  if (!is.null(Sigma_pop)) {
    sn <- intersect(colnames(Sigma_pop), names(pop_full))
    out$meat1 <- as.numeric(crossprod(pop_full[sn], Sigma_pop[sn, sn] %*% pop_full[sn]))
  }
  # target-score meat (level 2)
  if (is.null(Gs)) {
    Gs <- tryCatch(
      target_residual_scores_v2(dat, beta_target, tau, lambda_gamma, h = h,
                                nuisance_control = nuisance_control %||% list()),
      error = function(e) NULL
    )
  }
  proj_star <- if (!is.null(Gs) && ncol(Gs) == length(pop_full)) {
    as.numeric(Gs %*% pop_full)
  } else NULL
  if (!is.null(proj_star)) {
    out$meat2 <- mean((proj_star - mean(proj_star))^2)
    out$L_k <- -(1 / sqrt(n)) * sum(proj_star)
  }
  # fitted-score meat with pop direction (level 3)
  proj0 <- as.numeric(G0 %*% pop_full)
  out$meat3 <- mean((proj0 - mean(proj0))^2)
  # fitted-score meat with sample exact inverse row (level 4)
  omega_exact <- tryCatch(as.numeric(solve_linear_pd_v2(Hf, e, name = "exact row")),
                          error = function(e) NULL)
  if (!is.null(omega_exact)) {
    proj_ex <- as.numeric(G0 %*% omega_exact)
    out$meat4 <- mean((proj_ex - mean(proj_ex))^2)
  }
  out
}

# -----------------------------------------------------------------------------
# TRUE-SUPPORT finite-sample variance diagnostics (ROUND-2, item A1/8).
#
# All are DIAGNOSTICS ONLY for the low-dimensional (oracle) fit; none is used to
# promote a correction to primary inference or selected by coverage.
#
#   CR0  : the fitted unsmoothed-score sandwich (as currently reported).
#   CR3  : Mancl-DeRouen style: cluster leverage h_i = g_i'(G'G)^{-1}g_i and
#          meat inflated by (1 - h_i)^{-2} (diagonal one-borrowing analogue).
#   KC   : Kauermann-Carroll style: meat inflated by (1 - h_i)^{-1/2}.
# Each returns an SE = sqrt(omega' meat omega / n) on the per-cluster scale.
# -----------------------------------------------------------------------------
cluster_leverage_v2 <- function(cluster_scores) {
  G <- cluster_scores
  Gc <- sweep(G, 2, colMeans(G))
  S <- crossprod(Gc) / nrow(G)
  Sinv <- tryCatch(solve(S, tol = 1e-12), error = function(e) NULL)
  if (is.null(Sinv)) return(rep(1, nrow(G)))
  h <- diag(Gc %*% Sinv %*% t(Gc)) / nrow(G)
  pmax(h, 1e-10)
}

finite_sample_sandwich_diagnostics_v2 <- function(cluster_scores, omega) {
  G <- cluster_scores
  n <- nrow(G)
  omega <- as.numeric(omega)
  proj <- as.numeric(G %*% omega)
  c0 <- mean((proj - mean(proj))^2)
  h <- cluster_leverage_v2(G)
  cr3 <- mean(((proj - mean(proj))^2) / spow(1 - h, 2))
  kc <- mean(((proj - mean(proj))^2) / spow(1 - h, 0.5))
  list(
    se_cr0 = sqrt(max(c0, 0) / n),
    se_cr3 = sqrt(max(cr3, 0) / n),
    se_kc = sqrt(max(kc, 0) / n),
    leverage_max = max(h)
  )
}

spow <- function(x, e) {
  out <- x^e
  out[!is.finite(out)] <- 1
  out
}

# -----------------------------------------------------------------------------
# Delete-one-cluster jackknife for the TRUE-SUPPORT (oracle-support) estimator.
# Fits the unpenalised profile on the oracle support and computes the one-step
# debiased coordinate, for every delete-one-cluster subsample. Diagnostic only.
# -----------------------------------------------------------------------------
oracle_jackknife_v2 <- function(dat, active, target_coords, tau, h, lambda_gamma,
                                coords_fit_control = list()) {
  keep <- sort(unique(c(as.integer(active),
                        match(target_coords, colnames(dat$X)))))
  subX <- dat$X[, keep, drop = FALSE]
  target_names <- intersect(target_coords, colnames(subX))
  fit_once <- function(rows) {
    f <- tryCatch(
      fit_profile_lasso_v2(
        y = dat$y[rows], X = subX[rows, , drop = FALSE],
        Z = dat$Z[rows, , drop = FALSE], cluster_id = dat$cluster_id[rows],
        tau = tau, h = h, lambda_beta = 0, lambda_gamma = lambda_gamma,
        penalty_factor = rep(0, ncol(subX)),
        control = modifyList(list(max_iter = 100L, beta_tol = 1e-6, kkt_tol = 1e-5),
                             coords_fit_control)
      ),
      error = function(e) NULL
    )
    if (is.null(f)) return(NULL)
    Hf <- f$components$hessian
    Hf <- (Hf + t(Hf)) / 2
    out <- setNames(rep(NA_real_, length(target_names)), target_names)
    for (nm in target_names) {
      k <- match(nm, colnames(subX))
      e <- numeric(ncol(subX)); e[k] <- 1
      om <- tryCatch(as.numeric(solve(Hf, e)), error = function(e) NULL)
      if (is.null(om)) next
      out[nm] <- f$beta[k] - sum(om * f$components$score)
    }
    out
  }
  full <- fit_once(seq_len(length(dat$y)))
  if (is.null(full)) return(list(jack_se = setNames(rep(NA_real_, length(target_names)), target_names),
                                 full = NULL))
  rows_by_cluster <- split(seq_len(length(dat$y)), dat$cluster_id)
  leave_out <- lapply(rows_by_cluster, function(ix) setdiff(seq_len(length(dat$y)), ix))
  jack <- do.call(rbind, lapply(leave_out, fit_once))
  n <- length(rows_by_cluster)
  jack_se <- sqrt((n - 1) / n * apply((jack - sweep(jack, 2, colMeans(jack)))^2, 2, sum))
  list(jack_se = setNames(as.numeric(jack_se), colnames(jack)),
       full = full, n_clusters = n)
}
