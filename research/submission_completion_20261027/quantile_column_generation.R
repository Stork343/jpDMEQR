submission_quantile_fnb <- function(X, y, tau) {
  invisible(loadNamespace("quantreg"))
  n <- length(y)
  p <- ncol(X)
  workspace <- numeric(10L * n)
  workspace[seq_len(n)] <- 1 - tau
  native <- .Fortran("rqfnb", as.integer(n), as.integer(p), a = as.double(t(X)),
    c = as.double(-y), rhs = as.double((1 - tau) * colSums(X)), d = rep(1, n),
    u = rep(1, n), beta = as.double(.99995), eps = as.double(1e-9),
    wn = workspace, wp = double((p + 3L) * p), nit = integer(3L), info = integer(1L), PACKAGE = "quantreg")
  if (native$info != 0L) stop("The quantreg Frisch Newton solver returned a nonzero numerical code.")
  list(beta = -native$wp[seq_len(p)], dual = native$wn[2L * n + seq_len(n)] - (1 - tau), iterations = native$nit)
}

submission_quantile_column_fit <- function(design, response, tau, penalty) {
  X <- assert_numeric_matrix_v2(design, "design")
  y <- as.numeric(response)
  p <- ncol(X)
  N <- nrow(X)
  if (length(y) != N || any(!is.finite(y))) stop("Invalid quantile regression response.")
  if (length(penalty) == 1L) penalty <- rep(penalty, p)
  if (length(penalty) != p || any(!is.finite(penalty)) || any(penalty < 0)) stop("Invalid quantile penalties.")
  score_zero <- tau - (y <= quantile(y, tau, names = FALSE))
  ranking <- abs(as.numeric(crossprod(X, score_zero)))
  selected <- union(which(penalty == 0), head(order(ranking - penalty, decreasing = TRUE), min(16L, p)))
  records <- list()
  for (iteration in seq_len(p)) {
    local_penalty <- penalty[selected]
    positive <- which(local_penalty > 0)
    augmented <- X[, selected, drop = FALSE]
    augmented_y <- y
    if (length(positive)) {
      rows <- diag(local_penalty, nrow = length(selected))[positive, , drop = FALSE]
      augmented <- rbind(augmented, rows, -rows)
      augmented_y <- c(y, rep(0, 2L * length(positive)))
    }
    fitted <- submission_quantile_fnb(augmented, augmented_y, tau)
    dual <- fitted$dual[seq_len(N)]
    gradient <- as.numeric(crossprod(X, dual))
    omitted <- setdiff(seq_len(p), selected)
    excess <- abs(gradient) - penalty
    violations <- omitted[excess[omitted] > 1e-7 * pmax(1, penalty[omitted])]
    records[[iteration]] <- data.frame(iteration = iteration, columns = length(selected), violations = length(violations))
    if (!length(violations)) break
    selected <- union(selected, head(violations[order(excess[violations], decreasing = TRUE)], max(16L, length(selected))))
  }
  if (length(violations)) stop("The quantile column programme exhausted the full coordinate budget.")
  beta <- setNames(numeric(p), colnames(X))
  beta[selected] <- fitted$beta
  residual <- y - as.numeric(X %*% beta)
  objective <- sum(residual * (tau - (residual < 0))) + sum(penalty * abs(beta))
  dual_objective <- sum(y * dual)
  gap <- objective - dual_objective
  violation <- max(c(0, abs(gradient) - penalty))
  if (any(!is.finite(c(beta, dual, objective, gap, violation))) ||
      min(dual) < tau - 1 - 1e-8 || max(dual) > tau + 1e-8 ||
      violation > 1e-6 * max(1, max(penalty)) || abs(gap) > 1e-6 * max(1, objective))
    stop("The full quantile regression primal-dual certificate failed.")
  list(coefficients = beta, tau = tau, residuals = residual, nit = fitted$iterations,
    penalty = penalty, solver_eps = 1e-9, objective = objective, dual_objective = dual_objective,
    dual = dual, primal_dual_gap = gap, maximum_dual_violation = violation,
    column_records = do.call(rbind, records), solver = "quantreg_FNB_certified_column_generation")
}
