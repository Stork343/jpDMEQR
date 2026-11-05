library(osqp)

cq_qp_path <- function(features, tau, alpha, config) {
  lambda <- cq_lambda(features, tau, alpha, config)
  X <- features$X
  N <- nrow(X)
  p <- ncol(X)
  n <- features$n
  clinical <- which(features$penalty_factor == 0)
  penalized <- which(features$penalty_factor == 1)
  decomposition <- qr(X[, clinical, drop = FALSE], tol = 1e-10)
  stopifnot(decomposition$rank == length(clinical))
  projection <- qr.coef(decomposition, X[, penalized, drop = FALSE])
  X[, penalized] <- qr.resid(decomposition, X[, penalized, drop = FALSE])
  X[, clinical] <- qr.Q(decomposition) * sqrt(N)
  d <- if (alpha > 0) length(penalized) else 0L
  equality <- cbind(cq_sparse(X), Diagonal(N), -Diagonal(N), Matrix(0, N, d, sparse = TRUE))
  nonnegative <- cbind(Matrix(0, 2L * N, p, sparse = TRUE), Diagonal(2L * N),
    Matrix(0, 2L * N, d, sparse = TRUE))
  A <- rbind(equality, nonnegative)
  lower <- c(features$y, numeric(2L * N))
  upper <- c(features$y, rep(Inf, 2L * N))
  if (d > 0L) {
    select <- sparseMatrix(i = seq_len(d), j = penalized, x = 1, dims = c(d, p))
    positive <- cbind(select, Matrix(0, d, 2L * N, sparse = TRUE), -Diagonal(d))
    negative <- cbind(-select, Matrix(0, d, 2L * N, sparse = TRUE), -Diagonal(d))
    A <- rbind(A, positive, negative)
    lower <- c(lower, rep(-Inf, 2L * d))
    upper <- c(upper, numeric(2L * d))
  }
  P <- sparseMatrix(i = penalized, j = penalized, x = rep(n * lambda[1L] * (1 - alpha), length(penalized)),
    dims = c(p + 2L * N + d, p + 2L * N + d))
  linear <- c(numeric(p), tau * features$weight, (1 - tau) * features$weight,
    rep(n * lambda[1L] * alpha, d))
  solver <- osqp(P = P, q = linear, A = cq_sparse(A), l = lower, u = upper,
    pars = osqpSettings(eps_abs = config$path_tolerance, eps_rel = config$path_tolerance,
      max_iter = config$path_max_iterations,
      polish = TRUE, verbose = FALSE, warm_start = TRUE, check_termination = 25L))
  beta <- matrix(NA_real_, p, length(lambda), dimnames = list(colnames(X), NULL))
  checks <- list()
  for (index in seq_along(lambda)) {
    linear <- c(numeric(p), tau * features$weight, (1 - tau) * features$weight,
      rep(n * lambda[index] * alpha, d))
    solver$Update(q = linear, Px = rep(n * lambda[index] * (1 - alpha), length(penalized)))
    fit <- solver$Solve()
    if (fit$info$status != "solved") stop(paste("OSQP", index, fit$info$status))
    coordinates <- fit$x[seq_len(p)]
    beta[, index] <- coordinates
    beta[clinical[decomposition$pivot], index] <- backsolve(qr.R(decomposition), coordinates[clinical] * sqrt(N))
    beta[clinical, index] <- beta[clinical, index] - drop(projection %*% beta[penalized, index])
    stopifnot(max(abs(drop(features$X %*% beta[, index]) - drop(X %*% coordinates))) < 1e-8)
    objective <- cq_objective(features, beta[, index], tau, lambda[index], alpha)
    checks[[index]] <- data.table(index, iterations = fit$info$iter,
      primal_residual = fit$info$pri_res, dual_residual = fit$info$dua_res,
      objective_error = abs(objective - fit$info$obj_val / n))
    stopifnot(abs(objective - fit$info$obj_val / n) < config$maximum_path_objective_difference)
  }
  list(beta = beta, lambda = lambda, alpha = alpha,
    iterations = vapply(checks, function(value) value$iterations, numeric(1)),
    diagnostics = rbindlist(checks), solver = "osqp_exact_path")
}
