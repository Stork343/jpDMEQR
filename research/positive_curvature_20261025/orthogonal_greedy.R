pc_orthogonal_path <- function(prepared, tau, sizes, gap_tolerance = 1e-8, max_iterations = 30000L,
                                solver = joint_profile_pdal_cpp) {
  p <- prepared$p
  N <- nrow(prepared$Xs)
  square_weight <- sqrt(prepared$inverse_size / prepared$n)
  weighted <- prepared$Xs * square_weight
  selected <- which(prepared$penalty_factor == 0)
  excluded <- which(prepared$penalty_factor > 0)
  gamma <- prepared$gamma0
  dual <- prepared$dual0
  fitted <- NULL
  if (length(selected)) {
    decomposition <- qr(weighted[, selected, drop = FALSE], tol = 1e-10)
    if (decomposition$rank != length(selected)) stop("Always included columns are linearly dependent.")
    selected <- selected[decomposition$pivot]
    basis <- qr.Q(decomposition)
    triangular <- qr.R(decomposition)
    basis_coefficient <- drop(triangular %*% prepared$theta0[selected])
    fitted <- solver(
      basis / square_weight, prepared$Z, prepared$y, prepared$row_cluster,
      prepared$inverse_size, prepared$Lambda, numeric(length(selected)), tau,
      basis_coefficient, gamma, dual, max_iterations, gap_tolerance
    )
    basis_coefficient <- fitted$theta
    gamma <- fitted$gamma
    dual <- fitted$dual
    projections <- crossprod(basis, weighted)
  } else {
    basis <- matrix(0, N, 0)
    triangular <- matrix(0, 0, 0)
    projections <- matrix(0, 0, p)
    basis_coefficient <- numeric()
  }
  theta <- numeric(p)
  if (length(selected)) theta[selected] <- drop(backsolve(triangular, basis_coefficient))
  coefficients <- matrix(0, p, max(sizes) + 1L)
  coefficients[, 1] <- theta / prepared$scale
  selected_path <- vector("list", max(sizes) + 1L)
  selected_path[[1]] <- selected
  diagnostics <- vector("list", max(sizes))
  column_variance <- colSums(weighted^2)
  for (step in seq_len(max(sizes))) {
    weighted_dual <- dual / square_weight
    residual_dual <- weighted_dual - drop(basis %*% crossprod(basis, weighted_dual))
    scores <- drop(crossprod(weighted, residual_dual)) / prepared$n
    residual_variance <- column_variance - colSums(projections^2)
    eligible <- excluded[residual_variance[excluded] > 1e-10]
    if (!length(eligible)) stop("The selected design has no independent extension.")
    candidate <- eligible[which.max(abs(scores[eligible]) / sqrt(residual_variance[eligible]))]
    old_dimension <- length(selected)
    coefficient <- projections[, candidate]
    new_basis <- weighted[, candidate] - drop(basis %*% coefficient)
    correction <- drop(crossprod(basis, new_basis))
    coefficient <- coefficient + correction
    new_basis <- new_basis - drop(basis %*% correction)
    norm <- sqrt(sum(new_basis^2))
    if (!is.finite(norm) || norm <= 1e-10) stop("The proposed orthogonal extension is singular.")
    new_basis <- new_basis / norm
    triangular <- rbind(cbind(triangular, coefficient), c(rep(0, old_dimension), norm))
    basis <- cbind(basis, new_basis)
    projections <- rbind(projections, drop(crossprod(weighted, new_basis)))
    selected <- c(selected, candidate)
    excluded <- setdiff(excluded, candidate)
    fitted <- solver(
      basis / square_weight, prepared$Z, prepared$y,
      prepared$row_cluster, prepared$inverse_size, prepared$Lambda,
      numeric(length(selected)), tau, c(basis_coefficient, 0), gamma, dual,
      max_iterations, gap_tolerance
    )
    basis_coefficient <- fitted$theta
    gamma <- fitted$gamma
    dual <- fitted$dual
    theta <- numeric(p)
    theta[selected] <- drop(backsolve(triangular, basis_coefficient))
    coefficients[, step + 1L] <- theta / prepared$scale
    selected_path[[step + 1L]] <- selected
    diagnostics[[step]] <- list(selected_coordinate = candidate, gap = fitted$duality_gap,
                                objective = fitted$primal, iterations = fitted$iterations,
                                conditional_design_variance = norm^2,
                                extension_bound = prepared$n * scores[candidate]^2 / (2 * norm^2))
  }
  list(coefficients = coefficients, selected_path = selected_path, diagnostics = diagnostics,
       fitted = fitted, theta = theta, gamma = gamma, dual = dual)
}

pc_fit_orthogonal_profile <- function(y, X, Z, cluster_id, tau, coordinates,
                                     Lambda = diag(ncol(Z)), penalty_factor = rep(1, ncol(X)),
                                     h = length(unique(cluster_id))^(-0.3), confidence = 0.95) {
  validate_tau_h_v2(tau, h)
  if (!is.numeric(y) || any(!is.finite(y))) stop("The response must be finite and numeric.")
  if (is.null(colnames(X)) || anyNA(colnames(X)) || any(!nzchar(colnames(X))) || anyDuplicated(colnames(X))) {
    stop("X must have unique nonempty column names.")
  }
  if (any(!penalty_factor %in% c(0, 1))) stop("Orthogonal search uses zero or one inclusion indicators.")
  if (length(confidence) != 1L || !is.finite(confidence) || confidence <= 0 || confidence >= 1) stop("Invalid confidence level.")
  initial <- pc_fit_greedy_initial(y, X, Z, cluster_id, tau, Lambda, penalty_factor,
                                   path_solver = pc_orthogonal_path)
  initial$ridge_metric <- "cluster_weighted_prediction"
  initial$tuning$rule <- "cluster_profile_risk_orthogonal_size"
  fitted <- pc_fit_infer(y, X, Z, cluster_id, tau, coordinates, Lambda,
                         penalty_factor, h = h, initial_fit = initial)
  started <- proc.time()[[3]]
  inference <- pc_inference_from_saved(list(y = y, X = X, Z = Z, cluster_id = cluster_id),
                                       fitted, tau, Lambda, confidence = confidence)
  inference$timing <- c(fitted$timing[c("initial", "geometry", "inference")],
                        score_inversion = proc.time()[[3]] - started)
  inference$timing <- c(inference$timing, total = sum(inference$timing))
  inference$initial <- fitted$initial
  inference$directions <- fitted$directions
  inference$beta_center <- fitted$beta_center
  inference$h <- fitted$h
  inference$confidence <- confidence
  inference$method_version <- "orthogonal_profile_score_20261025"
  inference
}
