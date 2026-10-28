profile_bandwidth_pair_scores_v2 <- function(y, X, Z, cluster_id, beta, tau, h,
                                             lambda_gamma = 1,
                                             nuisance_control = list()) {
  at_h <- profile_components_v2(
    y, X, Z, cluster_id, beta, tau, h, lambda_gamma,
    need_hessian = FALSE, nuisance_control = nuisance_control,
    need_fold_sums = FALSE
  )
  at_half <- profile_components_v2(
    y, X, Z, cluster_id, beta, tau, h / 2, lambda_gamma,
    need_hessian = FALSE, nuisance_control = nuisance_control,
    need_fold_sums = FALSE
  )
  if (!isTRUE(at_h$nuisance_converged) || !isTRUE(at_half$nuisance_converged)) {
    stop("Both nuisance profiles must converge before computing target scores.")
  }
  2 * at_half$cluster_scores - at_h$cluster_scores
}

inferential_cluster_scores_v2 <- function(components) {
  if (identical(components$score_definition,
                 "profile_bandwidth_pair_2_half_minus_full")) {
    return(components$cluster_scores)
  }
  components$unsmoothed_cluster_scores
}

profile_bandwidth_pair_v2 <- function(y, X, Z, cluster_id, beta, tau, h,
                                      lambda_gamma = 1,
                                      nuisance_control = list(),
                                      identity_check = FALSE,
                                      assembly = NULL,
                                      need_fold_sums = TRUE) {
  validate_tau_h_v2(tau, h)
  X <- assert_numeric_matrix_v2(X, "X")
  if (is.null(colnames(X)) || anyDuplicated(colnames(X))) {
    stop("X must have unique column names.")
  }

  # 两个 bandwidth 分别求解完整的 nuisance profile。
  at_h <- profile_components_v2(
    y, X, Z, cluster_id, beta, tau, h, lambda_gamma = lambda_gamma,
    need_hessian = TRUE, identity_check = identity_check,
    nuisance_control = nuisance_control, assembly = assembly,
    need_fold_sums = need_fold_sums
  )
  at_half_h <- profile_components_v2(
    y, X, Z, cluster_id, beta, tau, h / 2, lambda_gamma = lambda_gamma,
    need_hessian = TRUE, identity_check = identity_check,
    nuisance_control = nuisance_control, assembly = assembly,
    need_fold_sums = need_fold_sums
  )
  if (!isTRUE(at_h$nuisance_converged) ||
      !isTRUE(at_half_h$nuisance_converged)) {
    stop("Every nuisance profile must converge at both bandwidths; gradients: ",
         signif(at_h$max_nuisance_gradient, 5), ", ",
         signif(at_half_h$max_nuisance_gradient, 5), ".")
  }
  if (!identical(dimnames(at_h$cluster_scores),
                 dimnames(at_half_h$cluster_scores))) {
    stop("Cluster and coordinate names differ between bandwidths.")
  }

  scores <- 2 * at_half_h$cluster_scores - at_h$cluster_scores
  hessian <- 2 * at_half_h$hessian - at_h$hessian
  fold_sums <- NULL
  fold_counts <- NULL
  if (isTRUE(need_fold_sums)) {
    if (!identical(at_h$hessian_fold_counts, at_half_h$hessian_fold_counts)) {
      stop("Cluster folds differ between bandwidths.")
    }
    fold_sums <- Map(function(a, b) 2 * a - b,
                     at_half_h$hessian_fold_sums, at_h$hessian_fold_sums)
    fold_counts <- at_h$hessian_fold_counts
  }
  identity <- if (isTRUE(identity_check)) {
    2 * at_half_h$hessian_identity - at_h$hessian_identity
  } else NULL

  list(
    objective = 2 * at_half_h$objective - at_h$objective,
    score = colMeans(scores),
    hessian = hessian,
    cluster_scores = scores,
    hessian_identity = identity,
    profile_identity_error = if (isTRUE(identity_check)) {
      max(abs(hessian - identity))
    } else NA_real_,
    hessian_fold_sums = fold_sums,
    hessian_fold_counts = fold_counts,
    n_clusters = nrow(scores),
    h = h,
    bandwidths = c(h / 2, h),
    bandwidth_weights = c(2, -1),
    score_definition = "profile_bandwidth_pair_2_half_minus_full",
    nuisance_converged = TRUE,
    max_nuisance_gradient = max(at_h$max_nuisance_gradient,
                                at_half_h$max_nuisance_gradient),
    hessian_min_eigenvalue = min(eigen(hessian, symmetric = TRUE,
                                      only.values = TRUE)$values),
    components_h = at_h,
    components_half_h = at_half_h
  )
}

profile_bandwidth_pair_wald_v2 <- function(components, beta_start, precision,
                                           coordinates, alpha = 0.05) {
  if (!identical(components$score_definition,
                 "profile_bandwidth_pair_2_half_minus_full")) {
    stop("A paired-bandwidth profile component object is required.")
  }
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a scalar in (0, 1).")
  }
  G <- assert_numeric_matrix_v2(components$cluster_scores, "cluster_scores")
  precision <- assert_numeric_matrix_v2(precision, "precision")
  beta_start <- as.numeric(beta_start)
  p <- ncol(G)
  if (nrow(G) < 2L || length(beta_start) != p ||
      any(!is.finite(beta_start)) || nrow(precision) != p) {
    stop("Finite beta_start and precision must match the profile design.")
  }
  if (!is.numeric(coordinates) || length(coordinates) == 0L ||
      anyNA(coordinates) || any(!is.finite(coordinates)) ||
      any(coordinates < 1L | coordinates > p) ||
      any(coordinates != as.integer(coordinates)) ||
      anyDuplicated(coordinates) || length(coordinates) != ncol(precision)) {
    stop("coordinates must identify each precision column uniquely.")
  }
  coordinates <- as.integer(coordinates)
  projected <- G %*% precision
  centred <- sweep(projected, 2, colMeans(projected))
  variance <- colMeans(centred^2)
  if (any(!is.finite(variance)) || any(variance <= 0)) {
    stop("Every projected cluster-score variance must be positive.")
  }
  estimate <- beta_start[coordinates] -
    as.numeric(crossprod(precision, colMeans(G)))
  se <- sqrt(variance / nrow(G))
  canonical <- diag(p)[, coordinates, drop = FALSE]
  inverse_defect <- apply(abs(components$hessian %*% precision - canonical),
                          2, max)
  critical <- stats::qnorm(1 - alpha / 2)
  data.frame(
    coordinate = colnames(G)[coordinates],
    beta_start = beta_start[coordinates],
    beta_tilde = estimate,
    se = se,
    ci_lower = estimate - critical * se,
    ci_upper = estimate + critical * se,
    inverse_defect = inverse_defect,
    n_clusters = nrow(G),
    score_definition = components$score_definition,
    stringsAsFactors = FALSE
  )
}
