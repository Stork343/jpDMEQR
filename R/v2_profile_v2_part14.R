profile_pair_precision_backend_v2 <- function(solver_opts = list()) {
  requested <- solver_opts$backend %||% Sys.getenv(
    "JPDMEQR_DANTZIG_SOLVER", unset = "clarabel_native"
  )
  if (requested %in% c("cvxr", "native", "clarabel_native")) return("clarabel_native")
  if (identical(requested, "dual_warm")) return("dual_warm")
  stop("Unsupported paired precision backend: ", requested)
}

profile_pair_solve_grid_v2 <- function(H, coordinate, mu_grid, solver_opts = list(),
                                       fit = NULL, context = list()) {
  if (any(!is.finite(mu_grid)) || any(mu_grid <= 0 | mu_grid >= 1)) {
    stop("Paired coordinate inference requires a grid within (0,1).")
  }
  options <- profile_pair_solver_options_v2(solver_opts)
  if (!identical(options$backend, "dual_warm")) {
    return(lapply(mu_grid, function(mu) profile_pair_solve_row_v2(
      H, coordinate, mu, solver_opts, fit
    )))
  }
  if (!requireNamespace("highs", quietly = TRUE)) {
    profile_pair_stop_v2("dependency", "The declared dual backend requires highs.", fit)
  }
  grid <- profile_pair_run_stage_v2(
    "precision_solver", fit,
    dantzig_dual_solve_grid_v2(H, coordinate, mu_grid,
                               solver_opts = list(control = options$highs_control))
  )
  accepted_status <- c("optimal", "unbounded", "unbounded_or_infeasible",
                        "infeasible_by_monotonicity")
  attempts <- data.frame(mu = grid$mu, status = grid$status,
                          residual = grid$residual, feasible = grid$feasible)
  if (any(!grid$status %in% accepted_status) ||
      any(grid$status == "optimal" & !grid$feasible)) {
    profile_pair_stop_v2(
      "precision_solver", "The dual Dantzig calculation failed its status or residual checks.",
      fit, c(context, list(coordinate = coordinate, attempts = attempts, grid = grid))
    )
  }
  p <- nrow(H)
  lapply(seq_along(mu_grid), function(j) {
    feasible <- isTRUE(grid$feasible[j])
    omega <- if (feasible) grid$omega[[j]] else rep(NA_real_, p)
    list(omega = omega, coordinate = as.integer(coordinate),
         coordinate_name = colnames(H)[coordinate] %||% as.character(coordinate),
         mu = mu_grid[j], residual = grid$residual[j],
         solver = "HiGHS-dual-simplex", solver_status = grid$status[j],
         l1_norm = if (feasible) sum(abs(omega)) else NA_real_,
         l2_norm = if (feasible) sqrt(sum(omega^2)) else NA_real_,
         attempts = attempts, feasible = feasible, solver_options = options)
  })
}
