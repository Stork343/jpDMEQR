soy_initial_basis_fit <- milk_reference_fit_small

soy_refined_basis_fit <- function(data, initial, tau, h, Lambda, weights = rep(1, data$n),
                                  ridge = data$n^(-2), gram = crossprod(data$X, data$X * data$weight) / data$n) {
  fitted <- soy_initial_basis_fit(data, initial, tau, h, Lambda, weights, ridge, gram)
  row_weight <- data$weight * weights[data$group]
  evaluate <- function(beta, gamma) {
    profile <- ds_profile(data, beta, tau, h, Lambda, gamma)
    list(profile = profile,
      gradient = -drop(crossprod(data$X, profile$psi * row_weight)) / data$n + ridge * drop(gram %*% beta),
      objective = sum(profile$objective * weights) / data$n + ridge * sum(beta * (gram %*% beta)) / 2)
  }
  beta <- fitted$beta
  current <- evaluate(beta, fitted$gamma)
  for (iteration in seq_len(20L)) {
    H <- sb_geometry_cpp(data$X, data$Z, current$profile$phi, data$order, data$pointers,
      Lambda, FALSE, weights) + ridge * gram
    factor <- chol(H)
    direction <- -drop(backsolve(factor, forwardsolve(t(factor), current$gradient)))
    if (max(abs(current$gradient)) <= 1e-11 && max(abs(direction)) <= 1e-9) break
    decrement <- -sum(current$gradient * direction)
    if (!is.finite(decrement) || decrement <= 0) stop("Invalid profile refinement direction.")
    accepted <- FALSE
    for (backtrack in 0:50) {
      step <- 2^(-backtrack)
      candidate <- evaluate(beta + step * direction, current$profile$gamma)
      allowance <- 64 * .Machine$double.eps * max(1, abs(current$objective))
      if (candidate$objective <= current$objective - 1e-4 * step * decrement + allowance) {
        accepted <- TRUE
        break
      }
    }
    if (!accepted) stop("The profile precision refinement failed its objective check.")
    beta <- beta + step * direction
    current <- candidate
    if (iteration == 20L) stop("The profile precision refinement did not converge.")
  }
  target <- c(1, rep(0, ncol(data$X) - 1L))
  inverse <- drop(backsolve(factor, forwardsolve(t(factor), target)))
  scores <- -as.numeric(rowsum(drop(data$X %*% inverse) * current$profile$psi * data$weight,
    data$group, reorder = FALSE))
  list(beta = beta, projected_scores = scores, gamma = current$profile$gamma,
    iterations = fitted$iterations + iteration - 1L, kkt = max(abs(current$gradient)),
    objective = current$objective, hessian = H, inverse = inverse)
}
