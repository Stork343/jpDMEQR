fit_lowdim_profile_pair_v2 <- function(y, X, Z, cluster_id, tau, h,
                                       lambda_gamma = 1, control = list()) {
  X <- assert_numeric_matrix_v2(X, "X")
  fit <- fit_profile_pair_v2(
    y, X, Z, cluster_id, tau, h_inf = h,
    lambda_gamma = lambda_gamma, penalty_factor = rep(0, ncol(X)),
    target_coordinates = colnames(X), control = control
  )
  J <- fit$components$hessian
  eigenvalues <- eigen((J + t(J)) / 2, symmetric = TRUE, only.values = TRUE)$values
  condition_number <- max(eigenvalues) / min(eigenvalues)
  if (min(eigenvalues) <= 1e-8 || !is.finite(condition_number) ||
      condition_number >= 1e10) {
    profile_pair_stop_v2("profile_geometry", "Reduced paired curvature failed.", fit,
                        list(eigenvalues = eigenvalues,
                             condition_number = condition_number))
  }
  precision <- solve(J)
  interval <- profile_bandwidth_pair_wald_v2(
    fit$components, fit$beta_refit, precision, seq_len(ncol(X))
  )
  fit$beta_pair <- setNames(interval$beta_tilde, colnames(X))
  fit$lowdim_diagnostics <- list(
    minimum_eigenvalue = min(eigenvalues), condition_number = condition_number,
    inverse_defect = max(abs(J %*% precision - diag(ncol(X))))
  )
  list(fit = fit, table = interval, precision = precision)
}

full_subject_profile_pair_jackknife_v2 <- function(y, X, Z, cluster_id, tau, h,
                                                   lambda_gamma = 1,
                                                   control = list()) {
  clusters <- cluster_index_v2(cluster_id)
  n <- clusters$n
  if (n < 4L) stop("At least four subjects are required for a paired jackknife.")
  full <- fit_lowdim_profile_pair_v2(y, X, Z, cluster_id, tau, h,
                                    lambda_gamma, control)
  rows <- split(seq_along(y), as.character(cluster_id))
  estimates <- matrix(NA_real_, nrow = n, ncol = ncol(X),
                      dimnames = list(names(rows), colnames(X)))
  diagnostics <- vector("list", n)
  for (i in seq_along(rows)) {
    keep <- setdiff(seq_along(y), rows[[i]])
    leave <- fit_lowdim_profile_pair_v2(
      y[keep], X[keep, , drop = FALSE], Z[keep, , drop = FALSE],
      cluster_id[keep], tau, h, lambda_gamma, control
    )
    estimates[i, ] <- leave$fit$beta_pair
    diagnostics[[i]] <- data.frame(
      omitted_subject = names(rows)[i], n_clusters = n - 1L,
      converged = leave$fit$converged,
      refit_gradient = leave$fit$post_refit_gradient_max,
      nuisance_gradient = leave$fit$components$max_nuisance_gradient,
      profile_identity_error = leave$fit$components$profile_identity_error,
      hessian_min_eigenvalue = leave$fit$lowdim_diagnostics$minimum_eigenvalue,
      hessian_condition_number = leave$fit$lowdim_diagnostics$condition_number,
      inverse_defect = leave$fit$lowdim_diagnostics$inverse_defect,
      stringsAsFactors = FALSE
    )
  }
  if (any(!is.finite(estimates))) stop("Every subject must have a valid paired leave-one fit.")
  centered <- sweep(estimates, 2, colMeans(estimates))
  list(fit = full$fit, table = full$table, precision = full$precision,
       estimates = estimates, se = sqrt((n - 1) / n * colSums(centered^2)),
       diagnostics = do.call(rbind, diagnostics), n_clusters = n,
       bandwidth_rule = "fixed_full_sample_bandwidth_in_leave_one_fits")
}

lowdim_profile_pair_inference_v2 <- function(y, X, Z, cluster_id, tau, h,
                                             lambda_gamma = 1, ci_level = 0.95,
                                             control = list()) {
  if (length(ci_level) != 1L || !is.finite(ci_level) || ci_level <= 0 || ci_level >= 1) {
    stop("Invalid confidence level.")
  }
  jack <- full_subject_profile_pair_jackknife_v2(
    y, X, Z, cluster_id, tau, h, lambda_gamma, control
  )
  critical <- stats::qnorm((1 + ci_level) / 2)
  beta <- as.numeric(jack$fit$beta_pair)
  sandwich_se <- jack$table$se
  table <- data.frame(
    coordinate = colnames(X), beta_hat = beta,
    beta_refit = as.numeric(jack$fit$beta_refit),
    score_source = "bandwidth_pair",
    score_definition = jack$fit$score_definition,
    method_version = jack$fit$method_version,
    target_type = "unsmoothed_regularised_profile", h = h, h_half = h / 2,
    tau = tau, ci_level = ci_level, n_clusters = jack$n_clusters,
    sandwich_se = sandwich_se,
    sandwich_lower = beta - critical * sandwich_se,
    sandwich_upper = beta + critical * sandwich_se,
    jackknife_se = jack$se,
    jackknife_lower = beta - critical * jack$se,
    jackknife_upper = beta + critical * jack$se,
    jackknife_completed = nrow(jack$estimates), stringsAsFactors = FALSE
  )
  list(table = table, fit = jack$fit, jackknife = jack,
       precision = jack$precision)
}
