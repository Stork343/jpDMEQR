milk_reference_fit_small <- ds_fit_small

milk_fit_small <- function(data, initial, tau, h, Lambda, weights = rep(1, data$n),
                          ridge = data$n^(-2), gram = crossprod(data$X, data$X * data$weight) / data$n) {
  dimension <- ncol(data$X)
  square_weight <- sqrt(data$weight / data$n)
  decomposition <- qr(data$X * square_weight, tol = 1e-12)
  if (decomposition$rank != dimension) stop("The selected empirical design is rank deficient.")
  triangular <- qr.R(decomposition)
  permutation <- decomposition$pivot
  transformed <- data
  transformed$X <- qr.Q(decomposition) / square_weight
  actual_gram <- crossprod(data$X, data$X * data$weight) / data$n
  if (max(abs(actual_gram - gram)) > 1e-10 * max(1, max(abs(gram))))
    stop("The declared prediction ridge differs from the observed design metric.")
  basis_initial <- drop(triangular %*% initial[permutation])
  fitted <- milk_reference_fit_small(transformed, basis_initial, tau, h, Lambda,
    weights = weights, ridge = ridge, gram = diag(dimension))
  beta <- inverse <- numeric(dimension)
  beta[permutation] <- drop(backsolve(triangular, fitted$beta))
  target <- as.numeric(permutation == 1L)
  basis_target <- drop(forwardsolve(t(triangular), target))
  basis_inverse <- drop(solve(fitted$hessian, basis_target))
  inverse[permutation] <- drop(backsolve(triangular, basis_inverse))
  profile <- ds_profile(transformed, fitted$beta, tau, h, Lambda, fitted$gamma)
  projection <- drop(transformed$X %*% basis_inverse)
  scores <- -as.numeric(rowsum(projection * profile$psi * data$weight,
    data$group, reorder = FALSE))
  row_weight <- data$weight * weights[data$group]
  gradient <- -drop(crossprod(data$X, profile$psi * row_weight)) / data$n + ridge * drop(gram %*% beta)
  prediction_error <- max(abs(drop(data$X %*% beta) - drop(transformed$X %*% fitted$beta)))
  if (max(abs(gradient)) > 1e-8 || prediction_error > 1e-7)
    stop("The original-coordinate prediction or gradient check failed.")
  hessian <- matrix(0, dimension, dimension)
  hessian[permutation, permutation] <- crossprod(triangular, fitted$hessian %*% triangular)
  list(beta = beta, projected_scores = scores, gamma = fitted$gamma,
    iterations = fitted$iterations, kkt = max(abs(gradient)), objective = fitted$objective,
    hessian = hessian, inverse = inverse, basis_hessian = fitted$hessian,
    prediction_error = prediction_error, basis_kkt = fitted$kkt)
}
