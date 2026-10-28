submission_clime_row <- function(H, coordinate, gamma, control, anchor, full_solver) {
  p <- nrow(H)
  target <- numeric(p)
  target[coordinate] <- 1
  largest <- which.max(abs(H[coordinate, ]))
  analytic <- numeric(p)
  analytic[largest] <- (1 - gamma) / H[coordinate, largest]
  analytic_residual <- max(abs(drop(H %*% analytic) - target))
  if (is.finite(analytic_residual) && analytic_residual <= gamma + 1e-9) {
    dual <- numeric(p)
    dual[coordinate] <- 1 / abs(H[coordinate, largest])
    return(list(feasible = TRUE, omega = analytic, residual = analytic_residual,
      objective = sum(abs(analytic)), dual_objective = sum(dual * target) - gamma * sum(abs(dual)),
      dual_residual = max(abs(drop(crossprod(H, dual)))) - 1,
      status = "optimal", certificate = "single_coordinate_primal_dual", models = 0L,
      history = data.frame(), support = largest))
  }
  if (is.null(anchor)) {
    solver <- full_solver()
    solver$clear_basis()
    cost <- dantzig_dual_cost_v2(p, coordinate, gamma)
    solver$L(seq_len(2L * p), cost)
    solver$solve(log_to_console = FALSE)
    status <- highs_model_status_label_v2(solver$status())
    if (status %in% c("unbounded", "unbounded_or_infeasible")) return(list(feasible = FALSE,
      omega = rep(NA_real_, p), residual = NA_real_, status = status,
      certificate = "complete_dual_programme", models = 1L, support = seq_len(p), info = solver$info()))
    if (status != "optimal") stop(paste("Complete CLIME programme status", status))
    solution <- solver$solution()
    omega <- dantzig_dual_omega_v2(solution$row_dual, p)
    residual <- max(abs(drop(H %*% omega) - target))
    dual_value <- -sum(cost * solution$col_value)
    if (residual > gamma + 1e-6 || abs(sum(abs(omega)) - dual_value) > 1e-6 * max(1, sum(abs(omega)))) stop("The complete CLIME programme certificate failed.")
    return(list(feasible = TRUE, omega = omega, residual = residual, status = status,
      objective = sum(abs(omega)), dual_objective = dual_value, certificate = "complete_dual_programme",
      models = 1L, support = seq_len(p), info = solver$info()))
  }
  anchor_cost <- sum(abs(anchor))
  anchor_column <- drop(H %*% anchor)
  if (!is.finite(anchor_cost) || anchor_cost <= 0 || max(abs(anchor_column - target)) > gamma + 1e-8) stop("The supplied CLIME anchor is not feasible.")
  selected <- sort(unique(c(coordinate, head(order(abs(H[coordinate, ]), decreasing = TRUE), 5L))))
  history <- list()
  repeat {
    d <- length(selected)
    reduced <- H[, selected, drop = FALSE]
    A <- cbind(reduced, -reduced, anchor_column, -anchor_column)
    cost <- c(rep(1, 2L * d), anchor_cost, anchor_cost)
    model <- highs::highs_model(L = cost, lower = rep(0, length(cost)), upper = rep(Inf, length(cost)),
      A = Matrix::Matrix(A, sparse = TRUE), lhs = target - gamma, rhs = target + gamma, maximum = FALSE)
    solver <- highs::highs_solver(model, control = control)
    solver$solve(log_to_console = FALSE)
    status <- highs_model_status_label_v2(solver$status())
    if (status != "optimal") stop(paste("Anchored CLIME programme status", status))
    solution <- solver$solution()
    omega <- anchor * (solution$col_value[2L * d + 1L] - solution$col_value[2L * d + 2L])
    omega[selected] <- omega[selected] + solution$col_value[seq_len(d)] - solution$col_value[d + seq_len(d)]
    dual <- solution$row_dual
    value <- sum(cost * solution$col_value)
    pricing <- abs(drop(crossprod(H, dual)))
    history[[length(history) + 1L]] <- data.frame(columns = d, objective = value,
      maximum_pricing = max(pricing), anchor_coefficient = solution$col_value[2L * d + 1L] - solution$col_value[2L * d + 2L])
    remaining <- setdiff(seq_len(p), selected)
    violating <- remaining[pricing[remaining] > 1 + 1e-8]
    if (!length(violating)) break
    selected <- sort(union(selected, head(violating[order(pricing[violating], decreasing = TRUE)], 25L)))
  }
  residual <- max(abs(drop(H %*% omega) - target))
  dual_value <- sum(dual * target) - gamma * sum(abs(dual))
  if (residual > gamma + 1e-6 || max(pricing) > 1 + 1e-7 ||
      abs(sum(abs(omega)) - dual_value) > 1e-6 * max(1, sum(abs(omega)))) stop("The full CLIME pricing certificate failed.")
  list(feasible = TRUE, status = "optimal", omega = omega, residual = residual,
    objective = sum(abs(omega)), dual_objective = dual_value, dual_residual = max(pricing) - 1,
    certificate = "complete_column_pricing", models = length(history), history = do.call(rbind, history), support = selected)
}

submission_clime_matrix_grid <- function(H, gamma_grid, solver_opts) {
  p <- nrow(H)
  gamma_grid <- sort(unique(gamma_grid))
  control <- highs::highs_control(threads = 1L, solver = "simplex", presolve = "on", log_to_console = FALSE,
    primal_feasibility_tolerance = 1e-9, dual_feasibility_tolerance = 1e-9)
  complete_solver <- NULL
  full_solver <- function() {
    if (is.null(complete_solver)) {
      model <- highs::highs_model(L = dantzig_dual_cost_v2(p, 1L, max(gamma_grid)),
        lower = rep(0, 2L * p), upper = rep(Inf, 2L * p), A = dantzig_dual_matrix_build_v2(H),
        lhs = rep(-Inf, 2L * p), rhs = rep(1, 2L * p), maximum = FALSE)
      complete_solver <<- highs::highs_solver(model, control = highs::highs_control(
        threads = 1L, solver = "simplex", presolve = "off", log_to_console = FALSE,
        simplex_strategy = 1L, simplex_scale_strategy = 2L,
        primal_feasibility_tolerance = 1e-8, dual_feasibility_tolerance = 1e-8))
    }
    complete_solver
  }
  grids <- lapply(seq_len(p), function(k) list(mu = gamma_grid, omega = vector("list", length(gamma_grid)),
    solver_records = vector("list", length(gamma_grid)), residual = rep(NA_real_, length(gamma_grid)),
    l1_norm = rep(NA_real_, length(gamma_grid)), status = rep("not_required_after_global_infeasibility", length(gamma_grid)),
    feasible = rep(NA, length(gamma_grid)), solver = "HiGHS-certified-column-generation", coordinate = k))
  records <- list()
  cutoff <- NA_real_
  models <- 0L
  for (index in order(gamma_grid, decreasing = TRUE)) {
    anchors <- submission_clime_anchors_cpp(H, gamma_grid[index])
    for (k in seq_len(p)) {
      anchor <- if (anchors$valid[k]) anchors$anchors[, k] else NULL
      fitted <- submission_clime_row(H, k, gamma_grid[index], control, anchor, full_solver)
      models <- models + fitted$models
      grids[[k]]$omega[[index]] <- fitted$omega
      grids[[k]]$solver_records[[index]] <- fitted
      grids[[k]]$residual[index] <- fitted$residual
      grids[[k]]$l1_norm[index] <- if (fitted$feasible) sum(abs(fitted$omega)) else NA_real_
      grids[[k]]$status[index] <- fitted$status
      grids[[k]]$feasible[index] <- fitted$feasible
      records[[length(records) + 1L]] <- data.frame(gamma = gamma_grid[index], coordinate = k,
        status = fitted$status, residual = fitted$residual, feasible = fitted$feasible,
        restricted_programmes = fitted$models, generated_columns = length(fitted$support), certificate = fitted$certificate)
      if (!fitted$feasible) {
        cutoff <- gamma_grid[index]
        break
      }
    }
    if (is.finite(cutoff)) break
  }
  list(grids = grids, attempts = do.call(rbind, records), global_infeasible_at = cutoff,
    model_constructions = models, solver_control = control, basis_policy = "certified_column_generation",
    feasibility_certificate = list(verified = FALSE, definition = "per_row_primal_dual_certificates"),
    source_precision_definition = "common_gamma_full_matrix_CLIME")
}
