soy_pooled_lp_path <- function(data, tau, nlambda = 40L, lambda_min = .01, lambda = NULL, progress = FALSE) {
  N <- length(data$y)
  p <- ncol(data$X)
  penalised <- which(data$penalty_factor > 0)
  fixed <- which(data$penalty_factor == 0)
  means <- numeric(p)
  means[penalised] <- colMeans(data$X[, penalised, drop = FALSE])
  X <- sweep(data$X, 2L, means, "-")
  scale <- sqrt(colMeans(X^2))
  if (any(!is.finite(scale)) || any(scale <= 1e-10)) stop("A pooled training column has zero scale.")
  X <- sweep(X, 2L, scale, "/")
  control <- highs::highs_control(threads = 1L, solver = "simplex", presolve = "off",
    simplex_strategy = 1L, simplex_scale_strategy = 2L, log_to_console = FALSE,
    primal_feasibility_tolerance = 1e-9, dual_feasibility_tolerance = 1e-9)
  make_solver <- function(selected, penalty, basis = NULL) {
    model <- highs::highs_model(L = data$y, lower = rep(tau - 1, N), upper = rep(tau, N),
      A = Matrix::Matrix(t(X[, selected, drop = FALSE]), sparse = TRUE),
      lhs = -penalty[selected], rhs = penalty[selected], maximum = TRUE)
    solver <- highs::highs_solver(model, control)
    if (!is.null(basis)) solver$set_basis(basis$col_status,
      c(basis$row_status, rep(1L, length(selected) - length(basis$row_status))))
    solver
  }
  selected <- fixed
  penalty <- numeric(p)
  solver <- make_solver(selected, penalty)
  solver$solve(log_to_console = FALSE)
  if (solver$status_message() != "Optimal") stop("The pooled unpenalised initial fit is not optimal: ", solver$status_message())
  initial <- solver$solution()
  gradient <- drop(crossprod(X, initial$col_value))
  if (is.null(lambda)) {
    maximum <- max(abs(gradient[penalised])) / N
    if (!is.finite(maximum) || maximum <= 0) stop("The pooled penalty entry value is not positive.")
    lambda <- maximum * exp(seq(0, log(lambda_min), length.out = nlambda))
  }
  if (any(!is.finite(lambda)) || any(lambda < 0) || any(diff(lambda) > 0)) stop("Invalid pooled penalty path.")
  coefficients <- matrix(0, p, length(lambda), dimnames = list(colnames(data$X), NULL))
  records <- list()
  for (point in seq_along(lambda)) {
    penalty <- N * lambda[point] * data$penalty_factor
    solver$cbounds(seq_along(selected), -penalty[selected], penalty[selected])
    generation <- 0L
    repeat {
      generation <- generation + 1L
      if (generation > p) stop("The pooled constraint generation exceeded its coordinate budget.")
      solver$solve(log_to_console = FALSE)
      if (solver$status_message() != "Optimal") stop("The pooled LP is not optimal: ", solver$status_message())
      solution <- solver$solution()
      dual <- solution$col_value
      gradient <- drop(crossprod(X, dual))
      remaining <- setdiff(seq_len(p), selected)
      excess <- abs(gradient) - penalty
      violating <- remaining[excess[remaining] / N > 1e-9]
      if (!length(violating)) break
      additional <- head(violating[order(excess[violating], decreasing = TRUE)], 50L)
      basis <- solver$get_basis()
      selected <- c(selected, additional)
      solver <- make_solver(selected, penalty, basis)
    }
    beta <- numeric(p)
    beta[selected] <- solution$row_dual
    residual <- data$y - drop(X %*% beta)
    primal <- sum(residual * (tau - (residual < 0))) + sum(penalty * abs(beta))
    dual_value <- sum(data$y * dual)
    gap <- (primal - dual_value) / N
    violation <- max(c(0, excess)) / N
    if (!is.finite(gap) || abs(gap) > 1e-7 || violation > 1e-8 ||
        min(dual) < tau - 1 - 1e-8 || max(dual) > tau + 1e-8)
      stop("The pooled path fails its complete primal-dual certificate.")
    original <- beta / scale
    intercept <- -sum(means * original)
    blocks <- grep("^Block[0-9]+$", colnames(data$X))
    original[blocks] <- original[blocks] + intercept
    identity <- max(abs(drop(data$X %*% original) - drop(X %*% beta)))
    if (identity > 1e-9) stop("The pooled scale transformation changed its predictions.")
    coefficients[, point] <- original
    records[[point]] <- data.table(point, lambda = lambda[point], active_constraints = length(selected),
      selected_variables = sum(abs(beta) > 1e-10), generations = generation,
      primal_dual_gap = gap, dual_violation = violation, mapping_error = identity)
    if (progress) {
      cat(format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE), "pooled penalty", point,
        "of", length(lambda), "constraints", length(selected), "gap", format(gap, digits = 3), "\n")
      flush.console()
    }
    if (point < length(lambda)) {
      basis <- solver$get_basis()
      keep <- selected %in% fixed | basis$row_status != 1L
      if (any(!keep)) {
        selected <- selected[keep]
        basis$row_status <- basis$row_status[keep]
        solver <- make_solver(selected, penalty, basis)
      }
    }
  }
  list(beta = coefficients, lambda = lambda, diagnostics = rbindlist(records),
    solver = "HiGHS_certified_pooled_quantile_path")
}
