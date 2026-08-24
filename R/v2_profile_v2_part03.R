fit_profile_lasso_v2 <- function(y,
                                 X,
                                 Z,
                                 cluster_id,
                                 tau,
                                 h,
                                 lambda_beta,
                                 lambda_gamma = 1,
                                 penalty_factor = NULL,
                                 beta_init = NULL,
                                 control = list()) {
  validate_tau_h_v2(tau, h)
  y <- as.numeric(y)
  X <- assert_numeric_matrix_v2(X, "X")
  Z <- assert_numeric_matrix_v2(Z, "Z")
  p <- ncol(X)
  if (!is.finite(lambda_beta) || lambda_beta < 0) stop("lambda_beta must be nonnegative.")

  penalty_factor <- penalty_factor %||% rep(1, p)
  penalty_factor <- as.numeric(penalty_factor)
  if (length(penalty_factor) != p || any(!is.finite(penalty_factor)) ||
      any(penalty_factor < 0)) stop("Invalid penalty_factor.")

  beta <- as.numeric(beta_init %||% rep(0, p))
  if (length(beta) != p) stop("beta_init has wrong length.")

  # Round-3 reference solver controls (METHOD_SPECIFICATION_ROUND3_AMENDMENT.md
  # section 3): max_iter >= 2000, max_backtrack >= 50, beta_tol = 1e-7.
  # The convergence beta-change rule is relative:
  #   last beta_change <= beta_tol * max(1, ||beta_hat||_inf),
  # and KKT acceptance uses the scale-normalised weighted residual R_KKT
  # (kkt_residual_normalized_v2) with R_KKT <= kkt_normalized_tol.
  max_iter <- control$max_iter %||% 2000L
  beta_tol <- control$beta_tol %||% 1e-7
  kkt_normalized_tol <- control$kkt_normalized_tol %||% 1e-3
  backtrack_factor <- control$backtrack_factor %||% 0.5
  max_backtrack <- control$max_backtrack %||% 50L
  min_step <- control$min_step %||% 1e-12
  verbose <- isTRUE(control$verbose)
  nuisance_control <- control$nuisance_control %||% list()

  effective_penalty <- lambda_beta * penalty_factor

  current <- profile_components_v2(
    y, X, Z, cluster_id, beta, tau, h, lambda_gamma,
    gamma_init = NULL, need_hessian = TRUE,
    nuisance_control = nuisance_control
  )
  composite <- current$objective + lambda_beta * sum(penalty_factor * abs(beta))
  history <- vector("list", max_iter)
  converged <- FALSE
  failure_stage <- ""

  # Optimiser implementation detail (statistical object unchanged: the
  # acceptance contract is the scale-normalised KKT). Accelerated proximal
  # gradient with momentum, monotone composite acceptance, and periodic
  # curvature (Hessian) refresh: intermediate candidates need only the
  # objective/score, which cuts the per-iteration cost by ~5x.
  accelerate <- isTRUE(control$accelerate %||% TRUE)
  refresh_every <- max(1L, as.integer(control$hessian_refresh %||% 5L))
  t_mom <- 1
  beta_prev <- beta
  step <- NA_real_

  for (iter in seq_len(max_iter)) {
    beta_prox <- if (accelerate && iter > 1L) {
      beta + ((t_mom - 1) / t_mom) * (beta - beta_prev)
    } else {
      beta
    }
    refresh <- !accelerate || (iter %% refresh_every) == 1L
    comp_prox <- tryCatch(
      profile_components_v2(
        y, X, Z, cluster_id, beta_prox, tau, h, lambda_gamma,
        gamma_init = current$gamma, need_hessian = refresh,
        nuisance_control = nuisance_control
      ),
      error = function(e) e
    )
    if (inherits(comp_prox, "error")) {
      failure_stage <- "profile_components"
      break
    }
    if (refresh) {
      eigmax <- largest_eigenvalue_power_v2(comp_prox$hessian)
      step <- 1 / max(eigmax, 1e-6)
    }
    score <- comp_prox$score
    accepted <- FALSE
    candidate <- NULL

    for (bt in 0:max_backtrack) {
      beta_new <- soft_threshold_weighted_v2(
        beta_prox - step * score,
        step * lambda_beta * penalty_factor
      )
      delta <- beta_new - beta_prox
      candidate <- tryCatch(
        profile_components_v2(
          y, X, Z, cluster_id, beta_new, tau, h, lambda_gamma,
          gamma_init = comp_prox$gamma, need_hessian = FALSE,
          nuisance_control = nuisance_control
        ),
        error = function(e) e
      )
      if (!inherits(candidate, "error")) {
        upper <- comp_prox$objective + sum(score * delta) +
          sum(delta^2) / (2 * step) + 1e-12
        composite_new <- candidate$objective +
          lambda_beta * sum(penalty_factor * abs(beta_new))
        if (candidate$objective <= upper && composite_new <= composite + 1e-12) {
          accepted <- TRUE
          break
        }
      }
      step <- step * backtrack_factor
      if (step < min_step) break
    }

    if (!accepted) {
      # Momentum restart: fall back to a plain proximal step from beta.
      if (accelerate && t_mom > 1) {
        t_mom <- 1
        beta_prev <- beta
        next
      }
      failure_stage <- "profile_backtracking"
      break
    }

    beta_change <- max(abs(beta_new - beta))
    kkt <- kkt_residual_profile_v2(
      beta_new, candidate$score, lambda_beta, penalty_factor
    )
    kkt_norm <- kkt_residual_normalized_v2(
      beta_new, candidate$score, lambda_beta * penalty_factor
    )

    history[[iter]] <- data.frame(
      iteration = iter,
      profile_objective = candidate$objective,
      composite_objective = composite_new,
      beta_change = beta_change,
      kkt_residual = kkt,
      kkt_normalized = kkt_norm,
      step = step,
      backtracks = bt,
      prox_refresh = refresh,
      nuisance_converged = candidate$nuisance_converged,
      max_nuisance_gradient = candidate$max_nuisance_gradient,
      profile_identity_error = candidate$profile_identity_error,
      stringsAsFactors = FALSE
    )

    if (verbose) {
      message(sprintf(
        "iter=%d composite=%.8f change=%.3e kkt=%.3e kkt_norm=%.3e step=%.3e",
        iter, composite_new, beta_change, kkt, kkt_norm, step
      ))
    }

    beta_prev <- beta
    beta <- beta_new
    current <- candidate
    composite <- composite_new
    if (accelerate) t_mom <- (1 + sqrt(1 + 4 * t_mom^2)) / 2

    beta_change_ok <- beta_change <= beta_tol * max(1, max(abs(beta_new)))
    if (beta_change_ok && kkt_norm <= kkt_normalized_tol &&
        current$nuisance_converged && all(is.finite(beta)) &&
        all(is.finite(candidate$score)) && is.finite(composite_new)) {
      converged <- TRUE
      break
    }
  }

  # The returned components must carry the Hessian (used by the precision and
  # variance stages); recompute once if the last accepted iterate was a
  # no-Hessian candidate.
  if (is.null(current$hessian)) {
    current <- tryCatch(
      profile_components_v2(
        y, X, Z, cluster_id, beta, tau, h, lambda_gamma,
        gamma_init = current$gamma, need_hessian = TRUE,
        nuisance_control = nuisance_control
      ),
      error = function(e) current
    )
  }

  history <- do.call(rbind, history[!vapply(history, is.null, logical(1))])
  final_kkt <- kkt_residual_profile_v2(
    beta, current$score, lambda_beta, penalty_factor
  )
  final_kkt_norm <- kkt_residual_normalized_v2(
    beta, current$score, lambda_beta * penalty_factor
  )

  list(
    beta = setNames(beta, colnames(X) %||% paste0("x", seq_len(p))),
    gamma = current$gamma,
    components = current,
    profile_objective = current$objective,
    composite_objective = composite,
    lambda_beta = lambda_beta,
    lambda_gamma = lambda_gamma,
    penalty_factor = penalty_factor,
    converged = converged,
    failure_stage = failure_stage,
    kkt_residual = final_kkt,
    kkt_normalized = final_kkt_norm,
    first_stage_beta_change = if (nrow(history)) tail(history$beta_change, 1) else NA_real_,
    first_stage_nonzero_count = sum(abs(beta) > 1e-10),
    iterations = nrow(history),
    history = history,
    call = match.call()
  )
}

# -----------------------------------------------------------------------------
# Round-3 cluster self-normalised first-stage calibration with round-4 dual
# bandwidths (METHOD_SPECIFICATION_ROUND4_AMENDMENT.md sections 1-5):
#   h_est = (log(p_P)/n)^{1/4}  -- first-stage estimation bandwidth only
#   h_inf = c_h n^{-3/10}       -- inferential bandwidth only
# All loading passes and first-stage fits run at h_est. After acceptance, the
# first-stage nuisance/score/Hessian objects are DISCARDED for inferential use;
# nuisance effects are reprofiled at h_inf and the exact inferential score,
# effective Hessian and fold Hessians are recomputed at the accepted beta_hat.
# The h_est components are retained ONLY in fit$components_first_stage for
# first-stage RSC diagnostics; they must never enter the Dantzig/one-step
# pipeline (guarded by assert_inferential_components_v2).
# -----------------------------------------------------------------------------
fit_profile_lasso_calibrated_v2 <- function(y,
                                            X,
                                            Z,
                                            cluster_id,
                                            tau,
                                            h_est = NULL,
                                            h_inf = NULL,
                                            lambda_gamma = 1,
                                            lambda_0_n,
                                            base_penalty_factor = NULL,
                                            beta_init = NULL,
                                            target_coordinates = NULL,
                                            control = list()) {
  X <- assert_numeric_matrix_v2(X, "X")
  p <- ncol(X)
  base_pf <- as.numeric(base_penalty_factor %||% rep(1, p))
  if (length(base_pf) != p || any(!is.finite(base_pf)) || any(base_pf < 0)) {
    stop("Invalid base_penalty_factor.")
  }
  P_idx <- which(base_pf > 0)
  n_pen <- length(P_idx)
  n_clusters_fit <- length(unique(cluster_id))
  # h_est fallback: (log(p_P)/n)^{1/4} from the FITTED design (round-4 1.1).
  if (is.null(h_est) || !is.finite(h_est) || h_est <= 0) {
    h_est <- (log(max(n_pen, 2)) / n_clusters_fit)^(1 / 4)
  }
  validate_tau_h_v2(tau, h_est)

  # No penalised coordinates: plain (unpenalised) fit, no loading machinery.
  # Such fits are inference-focused refits (e.g. TRUE-SUPPORT), so use h_inf.
  if (!n_pen) {
    fit <- fit_profile_lasso_v2(
      y, X, Z, cluster_id, tau, h_inf %||% h_est, lambda_beta = 0,
      lambda_gamma = lambda_gamma,
      penalty_factor = base_pf,
      beta_init = beta_init,
      control = modifyList(
        list(max_iter = 2000L, max_backtrack = 50L, beta_tol = 1e-7,
             kkt_normalized_tol = 1e-3),
        control$fit_control %||% list()
      )
    )
    fit$first_stage_calibrated <- FALSE
    return(fit)
  }
  if (is.null(h_inf) || !is.finite(h_inf) || h_inf <= 0) {
    stop("h_inf is required for the calibrated first stage (round-4 dual bandwidth).")
  }

  # lambda_0,n is a function of (n, p_P) of the FITTED problem; compute it from
  # the actual fit design (authoritative). For full-design methods this equals
  # tuning$lambda_0_n; screened sub-designs get their own value.
  lambda_meta <- lambda_0_n_v2(n_clusters_fit, n_pen)
  if (is.null(lambda_0_n)) {
    lambda_0_n <- lambda_meta$lambda_0
  } else if (abs(lambda_meta$lambda_0 - lambda_0_n) > 1e-12 * max(1, lambda_0_n)) {
    warning("lambda_0_n differs from the fit-design value; using the fit-design value.")
    lambda_0_n <- lambda_meta$lambda_0
  }
  if (!is.finite(lambda_0_n) || lambda_0_n <= 0) stop("lambda_0_n must be positive.")

  fit_control <- modifyList(
    list(max_iter = 2000L, max_backtrack = 50L, beta_tol = 1e-7,
         kkt_normalized_tol = 1e-3),
    control$fit_control %||% list()
  )
  nuisance_control <- control$nuisance_control %||%
    list(reltol = 1e-11, maxit = 400L, grad_tol = 1e-8)

  base_audit <- list(
    lambda_rule = "cluster-score-loading-v3",
    lambda_alpha = lambda_meta$alpha,
    lambda_normal_quantile = lambda_meta$normal_quantile,
    lambda_safety_constant = lambda_meta$safety_constant,
    lambda_base = lambda_0_n
  )

  loading_failure <- function(reason) {
    out <- base_audit
    out$beta <- setNames(rep(NA_real_, p), colnames(X))
    out$gamma <- NULL
    out$components <- NULL
    out$converged <- FALSE
    out$failure_stage <- reason
    out$first_stage_calibrated <- TRUE
    out$kkt_residual <- NA_real_
    out$kkt_normalized <- NA_real_
    out$iterations <- NA_integer_
    out
  }

  # Pass 0: loadings at b^(0) = 0 (nuisance profiled at 0), at h_est.
  comp0 <- tryCatch(
    profile_components_v2(
      y, X, Z, cluster_id, rep(0, p), tau, h_est, lambda_gamma,
      gamma_init = NULL, need_hessian = FALSE,
      nuisance_control = nuisance_control
    ),
    error = function(e) e
  )
  if (inherits(comp0, "error")) return(loading_failure("lambda_loading_pass0"))
  ell0 <- score_loadings_v2(comp0$cluster_scores)
  zero_profile_score_max <- max(abs(comp0$score))
  p0 <- lambda_0_n * ell0
  zero_kkt_ratio_max <- if (n_pen && any(p0[P_idx] > 0)) {
    max(abs(comp0$score[P_idx]) / p0[P_idx])
  } else NA_real_

  # Degenerate loadings (amendment 1.3): ell_ref = median positive loading on P.
  check_loadings <- function(ell) {
    ell_pos <- ell[P_idx]
    ell_ref <- stats::median(ell_pos[ell_pos > 0])
    !is.finite(ell_ref) || any(!is.finite(ell_pos)) ||
      any(ell_pos < 1e-8 * ell_ref)
  }
  if (check_loadings(ell0)) {
    return(loading_failure("lambda_loading_degenerate"))
  }

  # Preliminary fit with p_j^(0) = lambda_0,n * ell_j^(0) on P, at h_est.
  fit1 <- tryCatch(
    fit_profile_lasso_v2(
      y, X, Z, cluster_id, tau, h_est, lambda_beta = lambda_0_n,
      lambda_gamma = lambda_gamma,
      penalty_factor = ell0 * base_pf,
      beta_init = beta_init,
      control = fit_control
    ),
    error = function(e) e
  )
  if (inherits(fit1, "error")) {
    out <- loading_failure("penalised_fit")
    out$preliminary_error <- conditionMessage(fit1)
    return(out)
  }
  preliminary_kkt_normalized <- fit1$kkt_normalized %||% NA_real_

  # Pass 1: loadings at the preliminary fit b^(1), at h_est.
  comp1 <- tryCatch(
    profile_components_v2(
      y, X, Z, cluster_id, fit1$beta, tau, h_est, lambda_gamma,
      gamma_init = fit1$gamma, need_hessian = FALSE,
      nuisance_control = nuisance_control
    ),
    error = function(e) e
  )
  if (inherits(comp1, "error")) {
    out <- loading_failure("lambda_loading_pass1")
    out$preliminary_kkt_normalized <- preliminary_kkt_normalized
    return(out)
  }
  ell1 <- score_loadings_v2(comp1$cluster_scores)
  ell_final <- pmax(ell0, ell1)
  if (check_loadings(ell_final)) {
    out <- loading_failure("lambda_loading_degenerate")
    out$preliminary_kkt_normalized <- preliminary_kkt_normalized
    return(out)
  }

  # Final fit: warm start at b^(1), penalty p_j = lambda_0,n * ell_j_final,
  # at h_est.
  fit <- tryCatch(
    fit_profile_lasso_v2(
      y, X, Z, cluster_id, tau, h_est, lambda_beta = lambda_0_n,
      lambda_gamma = lambda_gamma,
      penalty_factor = ell_final * base_pf,
      beta_init = fit1$beta,
      control = fit_control
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    out <- loading_failure("penalised_fit")
    out$preliminary_kkt_normalized <- preliminary_kkt_normalized
    return(out)
  }

  # Acceptance contract (amendment section 3):
  #   R_KKT <= 1e-3; max_nuisance_gradient <= 1e-7;
  #   last beta_change <= 1e-7 * max(1, ||beta_hat||_inf); all finite.
  beta_change_ok <- is.finite(fit$first_stage_beta_change) &&
    fit$first_stage_beta_change <= 1e-7 * max(1, max(abs(fit$beta)))
  accepted <- isTRUE(fit$converged) &&
    is.finite(fit$kkt_normalized) && fit$kkt_normalized <= 1e-3 &&
    is.finite(fit$components$max_nuisance_gradient) &&
    fit$components$max_nuisance_gradient <= 1e-7 &&
    beta_change_ok && all(is.finite(fit$beta)) &&
    is.finite(fit$composite_objective)
  if (!accepted) {
    fit$converged <- FALSE
    fit$failure_stage <- "penalised_fit"
  }

  # Round-4 hard separation (amendment section 5): the h_est first-stage
  # components are diagnostics only; the inferential objects are recomputed at
  # h_inf at the accepted beta_hat.
  fit$components_first_stage <- fit$components  # h_est, diagnostics only

  # ---------------------------------------------------------------------------
  # Round-5 post-L1 profile refit (METHOD_SPECIFICATION_ROUND5_AMENDMENT.md
  # sections 1-5): the penalised L1 estimator is the SELECTION estimator only;
  # the inferential starting estimator is the zero-L1 profile refit on
  # S_R = S_hat U U U T at h_inf. Explicit failure; no fallback to beta_l1.
  # ---------------------------------------------------------------------------
  fit$beta_l1 <- fit$beta
  fit$post_refit_status <- "not_attempted"
  if (isTRUE(fit$converged)) {
    S_hat <- P_idx[abs(fit$beta[P_idx]) > 1e-10]
    T_idx <- if (length(target_coordinates) && !all(is.na(target_coordinates))) {
      sort(unique(match(as.character(target_coordinates), colnames(X))))
    } else integer()
    T_idx <- T_idx[!is.na(T_idx)]
    S_R <- sort(unique(c(S_hat, which(base_pf == 0), T_idx)))
    d_R <- length(S_R)
    refit_gate <- floor(n_clusters_fit / log(max(n_clusters_fit, 3)))
    fit$selected_support_size <- length(S_hat)
    fit$refit_set_size <- d_R
    fit$refit_contains_targets <- all(T_idx %in% S_R)
    if (d_R > refit_gate) {
      fit$converged <- FALSE
      fit$failure_stage <- "post_refit_dimension"
      fit$post_refit_status <- "dimension_gate_failed"
    } else {
      refit_control <- modifyList(
        list(max_iter = 1000L, max_backtrack = 50L, beta_tol = 1e-8,
             kkt_normalized_tol = 1e-7),
        control$refit_control %||% list()
      )
      refit <- tryCatch(
        fit_profile_lasso_v2(
          y = y, X = X[, S_R, drop = FALSE], Z = Z, cluster_id = cluster_id,
          tau = tau, h = h_inf, lambda_beta = 0, lambda_gamma = lambda_gamma,
          penalty_factor = rep(0, d_R), beta_init = fit$beta[S_R],
          control = refit_control
        ),
        error = function(e) e
      )
      if (inherits(refit, "error")) {
        fit$converged <- FALSE
        fit$failure_stage <- "post_refit"
        fit$post_refit_status <- "error"
      } else {
        H_R <- refit$components$hessian
        ev_R <- tryCatch(eigen((H_R + t(H_R)) / 2, symmetric = TRUE,
                               only.values = TRUE)$values,
                         error = function(e) numeric())
        grad_ok <- is.finite(refit$kkt_normalized) && refit$kkt_normalized <= 1e-7
        nuis_ok <- is.finite(refit$components$max_nuisance_gradient) &&
          refit$components$max_nuisance_gradient <= 1e-7
        bch_ok <- is.finite(refit$first_stage_beta_change) &&
          refit$first_stage_beta_change <= 1e-8 * max(1, max(abs(refit$beta)))
        eig_ok <- length(ev_R) && is.finite(min(ev_R)) && min(ev_R) > 1e-8 &&
          is.finite(max(ev_R) / min(ev_R)) && (max(ev_R) / min(ev_R)) < 1e10
        fin_ok <- all(is.finite(refit$beta)) && is.finite(refit$composite_objective)
        post_ok <- isTRUE(refit$converged) && grad_ok && nuis_ok && bch_ok &&
          eig_ok && fin_ok
        if (!post_ok) {
          fit$converged <- FALSE
          fit$failure_stage <- "post_refit"
          fit$post_refit_status <- "acceptance_failed"
        } else {
          beta_refit <- setNames(rep(0, p), colnames(X))
          beta_refit[S_R] <- refit$beta
          fit$beta_refit <- beta_refit
          fit$beta <- beta_refit
          fit$refit_set <- S_R
          fit$refit_set_names <- colnames(X)[S_R]
          fit$post_refit_status <- "ok"
          fit$post_refit_iterations <- refit$iterations
          fit$post_refit_gradient_max <- refit$kkt_normalized
          fit$post_refit_nuisance_gradient_max <- refit$components$max_nuisance_gradient
          fit$post_refit_beta_change <- refit$first_stage_beta_change
          fit$post_refit_hessian_min_eigenvalue <- min(ev_R)
          fit$post_refit_hessian_condition_number <- max(ev_R) / min(ev_R)
          # Full-p inferential reprofile at h_inf at beta_refit: every nuisance
          # profile, exact score, full effective Hessian and fold Hessians.
          # Reduced-model objects are discarded for the primary precision step.
          comp_final <- tryCatch(
            profile_components_v2(
              y, X, Z, cluster_id, beta_refit, tau, h_inf, lambda_gamma,
              gamma_init = NULL, need_hessian = TRUE,
              nuisance_control = nuisance_control
            ),
            error = function(e) e
          )
          if (!inherits(comp_final, "error")) {
            fit$components <- comp_final
            fit$gamma <- comp_final$gamma
          }
        }
      }
    }
  }

  p_final <- lambda_0_n * ell_final
  fit$first_stage_calibrated <- TRUE
  fit$lambda_rule <- base_audit$lambda_rule
  fit$lambda_alpha <- base_audit$lambda_alpha
  fit$lambda_normal_quantile <- base_audit$lambda_normal_quantile
  fit$lambda_safety_constant <- base_audit$lambda_safety_constant
  fit$lambda_base <- lambda_0_n
  fit$lambda_loading_pass0_min <- min(ell0[P_idx])
  fit$lambda_loading_pass0_median <- stats::median(ell0[P_idx])
  fit$lambda_loading_pass0_max <- max(ell0[P_idx])
  fit$lambda_loading_pass1_min <- min(ell1[P_idx])
  fit$lambda_loading_pass1_median <- stats::median(ell1[P_idx])
  fit$lambda_loading_pass1_max <- max(ell1[P_idx])
  fit$lambda_coordinate_min <- min(p_final[P_idx])
  fit$lambda_coordinate_median <- stats::median(p_final[P_idx])
  fit$lambda_coordinate_max <- max(p_final[P_idx])
  fit$zero_profile_score_max <- zero_profile_score_max
  fit$zero_kkt_ratio_max <- zero_kkt_ratio_max
  fit$preliminary_kkt_normalized <- preliminary_kkt_normalized
  fit$final_kkt_normalized <- fit$kkt_normalized
  fit$final_kkt_absolute <- fit$kkt_residual
  fit$first_stage_iterations <- fit$iterations
  fit$first_stage_beta_change <- fit$first_stage_beta_change
  fit$first_stage_nonzero_count <- fit$first_stage_nonzero_count
  fit$h_est <- h_est
  fit$h_inf <- h_inf
  # Round-5 post-refit audit fields (amendment section 10).
  fit$selected_support_size <- fit$selected_support_size %||% NA_integer_
  fit$refit_set_size <- fit$refit_set_size %||% NA_integer_
  fit$refit_contains_targets <- fit$refit_contains_targets %||% NA
  fit$post_refit_status <- fit$post_refit_status %||% "not_attempted"
  fit$post_refit_iterations <- fit$post_refit_iterations %||% NA_integer_
  fit$post_refit_gradient_max <- fit$post_refit_gradient_max %||% NA_real_
  fit$post_refit_nuisance_gradient_max <- fit$post_refit_nuisance_gradient_max %||% NA_real_
  fit$post_refit_beta_change <- fit$post_refit_beta_change %||% NA_real_
  fit$post_refit_hessian_min_eigenvalue <- fit$post_refit_hessian_min_eigenvalue %||% NA_real_
  fit$post_refit_hessian_condition_number <- fit$post_refit_hessian_condition_number %||% NA_real_
  fit
}

# Round-4 guard: the inferential precision/one-step pipeline must run on the
# h_inf-reprofiled components. This check fails if an h_est (or otherwise
# mismatched) Hessian/nuisance object is passed downstream.
assert_inferential_components_v2 <- function(fit, h_inf, tol = 1e-12) {
  if (isTRUE(fit$first_stage_calibrated)) {
    comp_h <- fit$components$h %||% NA_real_
    if (!is.finite(comp_h) || abs(comp_h - h_inf) > tol * max(1, h_inf)) {
      stop(sprintf(
        "Inferential components bandwidth mismatch: got h=%.6f, required h_inf=%.6f. ",
        comp_h, h_inf
      ))
    }
    if (is.null(fit$components_first_stage) ||
        !is.finite(fit$components_first_stage$h)) {
      stop("Calibrated fit is missing the h_est first-stage diagnostic components.")
    }
  }
  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# Row-wise Dantzig direction and debiased inference
# -----------------------------------------------------------------------------

solve_dantzig_row_v2 <- function(H,
                                 coordinate,
                                 mu_grid,
                                 solver_preference = c("CLARABEL", "ECOS", "SCS"),
                                 solver_opts = list()) {
  if (!requireNamespace("CVXR", quietly = TRUE)) {
    stop("Package 'CVXR' is required for the Dantzig program.")
  }
  H <- assert_numeric_matrix_v2(H, "H")
  H <- (H + t(H)) / 2
  p <- nrow(H)
  if (ncol(H) != p) stop("H must be square.")

  if (is.character(coordinate)) {
    if (is.null(colnames(H)) || !coordinate %in% colnames(H)) {
      stop("Unknown coordinate name: ", coordinate)
    }
    k <- match(coordinate, colnames(H))
  } else {
    k <- as.integer(coordinate)
    if (length(k) != 1L || k < 1L || k > p) stop("Invalid coordinate index.")
  }

  mu_grid <- sort(unique(as.numeric(mu_grid)))
  if (length(mu_grid) == 0L || any(!is.finite(mu_grid)) || any(mu_grid <= 0)) {
    stop("mu_grid must contain positive finite values.")
  }

  installed <- CVXR::installed_solvers()
  solver <- solver_preference[solver_preference %in% installed][1]
  if (is.na(solver)) stop("No requested CVXR solver is installed.")

  e <- numeric(p)
  e[k] <- 1
  attempts <- vector("list", length(mu_grid))

  for (ii in seq_along(mu_grid)) {
    mu <- mu_grid[ii]
    omega_var <- CVXR::Variable(p)
    constraints <- list(
      H %*% omega_var - e <= mu,
      H %*% omega_var - e >= -mu
    )
    problem <- CVXR::Problem(
      CVXR::Minimize(CVXR::p_norm(omega_var, 1)),
      constraints
    )
    # CVXR >= 1.9 renamed the exported solver entry point from solve to
    # psolve and changed the return contract: psolve returns the optimal
    # objective value, with the solution/status obtained via value() and
    # status(). Older CVXR returned a result object with $status/$getValue.
    cvxr_solve <- if ("psolve" %in% getNamespaceExports("CVXR")) {
      getExportedValue("CVXR", "psolve")
    } else {
      getExportedValue("CVXR", "solve")
    }
    solved <- tryCatch(
      do.call(cvxr_solve, c(list(problem = problem, solver = solver), solver_opts)),
      error = function(e) e
    )
    if (inherits(solved, "error")) {
      attempts[[ii]] <- data.frame(mu = mu, status = "solver_error",
                                   residual = NA_real_, message = conditionMessage(solved))
      next
    }
    if (is.numeric(solved) && length(solved) == 1L &&
        "value" %in% getNamespaceExports("CVXR") &&
        "status" %in% getNamespaceExports("CVXR")) {
      # CVXR >= 1.9 contract: solved is the objective value.
      status <- tryCatch(CVXR::status(problem), error = function(e) "unknown")
      omega <- tryCatch(as.numeric(CVXR::value(omega_var)), error = function(e) NULL)
    } else {
      status <- solved$status %||% "unknown"
      omega <- tryCatch(as.numeric(solved$getValue(omega_var)), error = function(e) NULL)
    }
    if (is.null(omega) || length(omega) != p || any(!is.finite(omega))) {
      attempts[[ii]] <- data.frame(mu = mu, status = status,
                                   residual = NA_real_, message = "No finite solution returned")
      next
    }
    residual <- max(abs(as.numeric(H %*% omega - e)))
    attempts[[ii]] <- data.frame(mu = mu, status = status,
                                 residual = residual, message = "")

    if (residual <= mu + 1e-6 && status %in% c("optimal", "optimal_inaccurate")) {
      return(list(
        omega = omega,
        coordinate = k,
        coordinate_name = colnames(H)[k] %||% as.character(k),
        mu = mu,
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
    mu = NA_real_,
    residual = NA_real_,
    solver = solver,
    solver_status = "infeasible_on_grid",
    l1_norm = NA_real_,
    l2_norm = NA_real_,
    attempts = do.call(rbind, attempts),
    feasible = FALSE
  )
}
