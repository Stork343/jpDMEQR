jp_fit_fixed <- function(y, X, Z, cluster_id, tau, coordinates = colnames(X),
                          Lambda = diag(ncol(Z)), h = length(unique(cluster_id))^(-0.3),
                          confidence = 0.95, initial = NULL) {
  validate_tau_h_v2(tau, h)
  X <- assert_numeric_matrix_v2(X, "X")
  Z <- assert_numeric_matrix_v2(Z, "Z")
  data <- sb_data(as.numeric(y), X, Z, cluster_id)
  Lambda <- make_lambda_matrix_v2(Lambda, ncol(Z))
  if (is.character(coordinates)) coordinates <- match(coordinates, colnames(X))
  if (!length(coordinates) || anyNA(coordinates) || anyDuplicated(coordinates) ||
      any(coordinates != as.integer(coordinates)) || any(coordinates < 1L | coordinates > ncol(X))) stop("Invalid fixed design coordinates.")
  square_weight <- sqrt(data$weight / data$n)
  decomposition <- qr(X * square_weight, tol = 1e-10)
  if (decomposition$rank != ncol(X)) stop("The fixed design is rank deficient.")
  if (is.null(initial)) initial <- drop(qr.coef(decomposition, y * square_weight))
  full <- ds_fit_small(data, initial, tau, h, Lambda)
  half <- ds_fit_small(data, full$beta, tau, h / 2, Lambda)
  target <- diag(ncol(X))[, coordinates, drop = FALSE]
  omega_full <- solve(full$hessian, target)
  omega_half <- solve(half$hessian, target)
  profile_full <- ds_profile(data, full$beta, tau, h, Lambda, full$gamma)
  profile_half <- ds_profile(data, half$beta, tau, h / 2, Lambda, half$gamma)
  Gfull <- -rowsum(X * (data$weight * as.numeric(profile_full$psi)), data$group, reorder = FALSE)
  Ghalf <- -rowsum(X * (data$weight * as.numeric(profile_half$psi)), data$group, reorder = FALSE)
  influence <- -2 * Ghalf %*% omega_half + Gfull %*% omega_full
  centered <- sweep(influence, 2L, colMeans(influence), "-")
  covariance <- crossprod(centered) / data$n^2
  se <- sqrt(diag(covariance))
  if (any(!is.finite(se)) || any(se <= 0)) stop("Nonpositive fixed design uncertainty.")
  beta <- setNames(2 * half$beta - full$beta, colnames(X))
  critical <- qnorm((1 + confidence) / 2)
  table <- data.frame(coordinate = colnames(X)[coordinates], beta_hat = beta[coordinates],
                       se = se, lower = beta[coordinates] - critical * se,
                       upper = beta[coordinates] + critical * se, confidence = confidence)
  diagnostics <- data.frame(
    n_clusters = data$n, n_visits = length(y), dimension = ncol(X), h = h,
    full_kkt = full$kkt, half_kkt = half$kkt,
    full_nuisance_gradient = max(profile_full$gradient_max),
    half_nuisance_gradient = max(profile_half$gradient_max),
    full_minimum_curvature = min(eigen(full$hessian, symmetric = TRUE, only.values = TRUE)$values),
    half_minimum_curvature = min(eigen(half$hessian, symmetric = TRUE, only.values = TRUE)$values),
    full_iterations = full$iterations, half_iterations = half$iterations)
  list(table = table, beta = beta, full = full, half = half, covariance = covariance,
       influence = influence, diagnostics = diagnostics,
       method_version = "joint_profile_fixed_design_20261026")
}
