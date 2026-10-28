submission_clime_primal_grid <- function(H, gamma_grid, solver_opts = list()) {
  p <- nrow(H)
  gamma_grid <- sort(unique(gamma_grid))
  control <- do.call(highs::highs_control, modifyList(list(
    threads = 1L, solver = "simplex", presolve = "off", log_to_console = FALSE,
    simplex_strategy = 1L, simplex_scale_strategy = 2L,
    primal_feasibility_tolerance = 1e-8, dual_feasibility_tolerance = 1e-8
  ), solver_opts$highs_options %||% list()))
  model <- highs::highs_model(L = rep(1, 2L * p), lower = rep(0, 2L * p),
    upper = rep(Inf, 2L * p), A = cbind(H, -H), lhs = rep(-1, p), rhs = rep(1, p))
  solver <- highs::highs_solver(model, control)
  parametric <- identical(solver_opts$algorithm, "parametric_warm")
  flare <- identical(solver_opts$algorithm, "flare_warm")
  grids <- lapply(seq_len(p), function(k) list(mu = gamma_grid,
    omega = vector("list", length(gamma_grid)), solver_records = vector("list", length(gamma_grid)),
    residual = rep(NA_real_, length(gamma_grid)), l1_norm = rep(NA_real_, length(gamma_grid)),
    status = rep("not_required_after_global_infeasibility", length(gamma_grid)),
    feasible = rep(NA, length(gamma_grid)), solver = "HiGHS-dual-simplex-primal-programme", coordinate = k))
  records <- list()
  cutoff <- NA_real_
  for (index in order(gamma_grid, decreasing = TRUE)) {
    gamma <- gamma_grid[index]
    warm_directions <- NULL
    for (k in seq_len(p)) {
      canonical <- numeric(p)
      canonical[k] <- 1
      # 目标函数保持不变，更新约束区间后继续使用满足dual条件的basis。
      solver$cbounds(seq_len(p), canonical - gamma, canonical + gamma)
      if (parametric && k > 1L) {
        warm <- submission_clime_parametric_basis(H, k, gamma)
        solver$set_basis(warm$col_status, warm$row_status)
      }
      if (flare && k > 1L) {
        if (is.null(warm_directions)) warm_directions <- submission_flare_matrix(H, gamma)
        omega_start <- warm_directions$matrices[[1L]][, k]
        omega_start[abs(omega_start) <= 10 * warm_directions$precision] <- 0
        warm <- submission_clime_basis(H, k, gamma, omega_start)
        solver$set_basis(warm$col_status, warm$row_status)
      }
      solver$solve(log_to_console = FALSE)
      status <- highs_model_status_label_v2(solver$status())
      grids[[k]]$status[index] <- status
      if (identical(status, "infeasible")) {
        grids[[k]]$feasible[index] <- FALSE
        cutoff <- gamma
        records[[length(records) + 1L]] <- data.frame(gamma = gamma, coordinate = k,
          status = status, residual = NA_real_, feasible = FALSE, dual_gap = NA_real_)
        break
      }
      if (!identical(status, "optimal")) {
        condition <- simpleError("The primal CLIME programme returned an unaccepted solver status.")
        condition$failure_stage <- "precision_solver"
        condition$diagnostics <- list(gamma = gamma, coordinate = k, status = status, info = solver$info())
        stop(condition)
      }
      solution <- solver$solution()
      omega <- solution$col_value[seq_len(p)] - solution$col_value[p + seq_len(p)]
      residual <- max(abs(H %*% omega - canonical))
      dual <- solution$row_dual
      dual_residual <- max(abs(crossprod(H, dual)))
      objective <- sum(abs(omega))
      dual_gap <- objective - dual[k] + gamma * sum(abs(dual))
      if (any(!is.finite(c(omega, residual, dual_gap, dual_residual))) || residual > gamma + 1e-6 ||
          dual_residual > 1 + 1e-6 || abs(dual_gap) > 1e-6 * max(1, objective)) {
        condition <- simpleError("The primal CLIME programme failed its primal-dual certificate.")
        condition$failure_stage <- "precision_solver"
        condition$diagnostics <- list(gamma = gamma, coordinate = k, residual = residual,
          dual_gap = dual_gap, dual_residual = dual_residual, solution = solution, info = solver$info())
        stop(condition)
      }
      grids[[k]]$omega[[index]] <- omega
      grids[[k]]$solver_records[[index]] <- list(info = solver$info(), dual_gap = dual_gap,
        dual_residual = dual_residual, dual = dual)
      grids[[k]]$residual[index] <- residual
      grids[[k]]$l1_norm[index] <- objective
      grids[[k]]$feasible[index] <- TRUE
      records[[length(records) + 1L]] <- data.frame(gamma = gamma, coordinate = k,
        status = status, residual = residual, feasible = TRUE, dual_gap = dual_gap)
    }
    if (is.finite(cutoff)) break
  }
  list(grids = grids, attempts = do.call(rbind, records), global_infeasible_at = cutoff,
    model_constructions = 1L, solver_control = control,
    basis_policy = if (flare) "flare_ADMM_and_HiGHS_primal_dual_certificate" else
      if (parametric) "fastclime_path_and_HiGHS_primal_dual_certificate" else "reuse_dual_feasible_with_row_bounds",
    feasibility_certificate = list(verified = FALSE, definition = "per_row_primal_dual_certificates"),
    source_precision_definition = "common_gamma_full_matrix_CLIME")
}
