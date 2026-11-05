cq_elastic_path <- function(features, tau, alpha, config) {
  stopifnot(alpha > 0, alpha < 1)
  lambda <- cq_lambda(features, tau, alpha, config)
  p <- ncol(features$X)
  clinical <- which(features$penalty_factor == 0)
  penalized <- which(features$penalty_factor == 1)
  beta <- matrix(0, p, length(lambda), dimnames = list(colnames(features$X), NULL))
  checks <- list()
  previous_gradient <- numeric(p)
  previous_beta <- numeric(p)
  for (index in seq_along(lambda)) {
    strong_threshold <- if (index == 1L) Inf else alpha * (2 * lambda[index] - lambda[index - 1L])
    selected <- sort(unique(c(clinical, which(abs(previous_beta) > 1e-9),
      which(abs(previous_gradient) > strong_threshold))))
    generation <- 0L
    iterations <- 0L
    repeat {
      generation <- generation + 1L
      stopifnot(generation <= ceiling(p / 64) + 2L)
      restricted <- features
      restricted$X <- features$X[, selected, drop = FALSE]
      restricted$penalty_factor <- features$penalty_factor[selected]
      fit <- cq_exact(restricted, tau, lambda[index], alpha, config)
      iterations <- iterations + fit$diagnostics$iterations
      current <- numeric(p)
      current[selected] <- fit$beta
      gradient <- drop(crossprod(features$X, fit$dual))
      violation <- abs(gradient) - lambda[index] * alpha
      omitted <- setdiff(penalized, selected)
      additions <- omitted[violation[omitted] > 1e-9]
      if (!length(additions)) break
      additions <- additions[order(violation[additions], decreasing = TRUE)]
      selected <- sort(unique(c(selected, head(additions, 64L))))
    }
    conjugate <- sum(pmax(abs(gradient[penalized]) - lambda[index] * alpha, 0)^2) /
      (2 * lambda[index] * (1 - alpha))
    objective <- cq_objective(features, current, tau, lambda[index], alpha)
    gap <- objective - sum(features$y * fit$dual) + conjugate
    stopifnot(abs(gap) < config$exact_fit_tolerance * max(1, abs(objective)))
    beta[, index] <- current
    previous_beta <- current
    previous_gradient <- gradient
    checks[[index]] <- data.table(index, iterations, generations = generation,
      included_columns = length(selected), gap,
      omitted_gradient_violation = if (length(omitted)) max(0, violation[omitted]) else 0)
  }
  list(beta = beta, lambda = lambda, alpha = alpha,
    iterations = vapply(checks, function(value) value$iterations, numeric(1)),
    diagnostics = rbindlist(checks), solver = "clarabel_certified_active_elastic_path")
}
