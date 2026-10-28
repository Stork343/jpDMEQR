submission_quantile_highs <- function(design, response, tau, penalty) {
  X <- assert_numeric_matrix_v2(design, "design")
  y <- as.numeric(response)
  N <- nrow(X)
  p <- ncol(X)
  if (length(y) != N || any(!is.finite(y))) stop("Invalid quantile response.")
  if (length(penalty) == 1L) penalty <- rep(penalty, p)
  if (length(penalty) != p || any(!is.finite(penalty)) || any(penalty < 0)) stop("Invalid quantile penalty.")
  model <- highs::highs_model(L = y, lower = rep(tau - 1, N), upper = rep(tau, N),
    A = t(X), lhs = -penalty, rhs = penalty, maximum = TRUE)
  solver <- highs::highs_solver(model, highs::highs_control(threads = 1L, solver = "simplex",
    simplex_strategy = 1L, simplex_scale_strategy = 2L, presolve = "on", log_to_console = FALSE,
    primal_feasibility_tolerance = 1e-8, dual_feasibility_tolerance = 1e-8))
  solver$solve(log_to_console = FALSE)
  status <- highs_model_status_label_v2(solver$status())
  if (status != "optimal") stop("The quantile dual LP did not return an optimal solution.")
  solution <- solver$solution()
  beta <- setNames(solution$row_dual, colnames(X))
  dual <- solution$col_value
  residual <- y - as.numeric(X %*% beta)
  objective <- sum(residual * (tau - (residual < 0))) + sum(penalty * abs(beta))
  dual_objective <- sum(y * dual)
  gap <- objective - dual_objective
  violation <- max(c(0, abs(as.numeric(crossprod(X, dual))) - penalty))
  if (any(!is.finite(c(beta, dual, gap))) ||
      violation > 1e-6 * max(1, max(penalty)) || abs(gap) > 1e-6 * max(1, objective))
    stop("The complete quantile dual LP failed its primal-dual certificate.")
  list(coefficients = beta, tau = tau, residuals = residual, penalty = penalty,
    objective = objective, dual_objective = dual_objective, dual = dual,
    primal_dual_gap = gap, maximum_dual_violation = violation, solver_eps = 1e-8,
    solver = "HiGHS_quantile_dual", solver_info = solver$info())
}
