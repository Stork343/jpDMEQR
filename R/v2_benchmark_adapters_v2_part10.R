iid_matrix_feasibility_certificate_v2 <- function(H, gamma_minimum) {
  eigenvalues <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  minimum <- min(eigenvalues)
  maximum <- max(eigenvalues)
  spectral_condition <- is.finite(minimum) && is.finite(maximum) &&
    maximum > 0 && minimum > 1e-12 * maximum
  certificate <- list(verified = FALSE, matrix_dimension = nrow(H),
                       gamma = gamma_minimum, minimum_eigenvalue = minimum,
                       maximum_eigenvalue = maximum, spectral_threshold = 1e-12,
                       spectral_condition = spectral_condition,
                       max_residual = NA_real_, coordinate_residuals = NULL,
                       direction_matrix = NULL, residual_matrix = NULL)
  if (!spectral_condition) return(certificate)
  inverse <- solve(H)
  residual <- profile_dense_crossprod_v2(H, inverse) - diag(nrow(H))
  certificate$direction_matrix <- inverse
  certificate$residual_matrix <- residual
  certificate$coordinate_residuals <- apply(abs(residual), 2L, max)
  certificate$max_residual <- max(certificate$coordinate_residuals)
  certificate$verified <- all(is.finite(inverse)) && is.finite(certificate$max_residual) &&
    certificate$max_residual <= gamma_minimum
  certificate
}

dantzig_iid_solve_matrix_grid_v2 <- function(H, gamma_grid, solver_opts = list()) {
  H <- assert_numeric_matrix_v2(H, "H")
  p <- nrow(H)
  gamma_grid <- sort(unique(as.numeric(gamma_grid)))
  if (ncol(H) != p || !length(gamma_grid) || any(!is.finite(gamma_grid)) ||
      any(gamma_grid <= 0 | gamma_grid >= 1)) stop("A square iid Hessian and tolerances in (0,1) are required.")
  H <- (H + t(H)) / 2
  if (identical(solver_opts$algorithm, "hybrid_certified")) {
    return(submission_clime_parametric_certified_grid(H, gamma_grid, modifyList(solver_opts, list(initializer = "hybrid"))))
  }
  if (identical(solver_opts$algorithm, "qp_certified")) {
    return(submission_clime_parametric_certified_grid(H, gamma_grid, modifyList(solver_opts, list(initializer = "qp"))))
  }
  if (identical(solver_opts$algorithm, "parametric_certified")) {
    return(submission_clime_parametric_certified_grid(H, gamma_grid, solver_opts))
  }
  if (identical(solver_opts$algorithm, "constraint_generation")) {
    return(submission_clime_constraint_grid(H, gamma_grid, solver_opts))
  }
  if (identical(solver_opts$algorithm, "certified_flare")) {
    return(submission_clime_certified_grid(H, gamma_grid, solver_opts))
  }
  if (identical(solver_opts$algorithm, "primal_warm")) {
    return(submission_clime_primal_grid(H, gamma_grid, solver_opts))
  }
  if (identical(solver_opts$algorithm, "parametric_warm")) {
    return(submission_clime_primal_grid(H, gamma_grid, solver_opts))
  }
  if (identical(solver_opts$algorithm, "flare_warm")) {
    return(submission_clime_primal_grid(H, gamma_grid, solver_opts))
  }
  if (identical(solver_opts$algorithm, "column_generation")) {
    return(submission_clime_matrix_grid(H, gamma_grid, solver_opts))
  }
  certificate <- iid_matrix_feasibility_certificate_v2(H, min(gamma_grid))
  primal_reuse <- identical(solver_opts$basis_policy, "reuse_primal_feasible")
  reset_basis <- !primal_reuse && !certificate$spectral_condition
  control <- do.call(highs::highs_control, modifyList(list(
    threads = 1L, solver = "simplex", presolve = "off", log_to_console = FALSE,
    simplex_strategy = if (primal_reuse) 4L else 1L,
    simplex_scale_strategy = if (reset_basis || primal_reuse) 2L else 0L,
    primal_feasibility_tolerance = 1e-8, dual_feasibility_tolerance = 1e-8
  ), solver_opts$highs_options %||% list()))
  model <- highs::highs_model(
    L = dantzig_dual_cost_v2(p, 1L, max(gamma_grid)),
    lower = rep(0, 2L * p), upper = rep(Inf, 2L * p),
    A = dantzig_dual_matrix_build_v2(H), lhs = rep(-Inf, 2L * p),
    rhs = rep(1, 2L * p), maximum = FALSE
  )
  solver <- highs::highs_solver(model, control = control)
  grids <- lapply(seq_len(p), function(k) list(
    mu = gamma_grid, omega = vector("list", length(gamma_grid)),
    solver_records = vector("list", length(gamma_grid)),
    residual = rep(NA_real_, length(gamma_grid)), l1_norm = rep(NA_real_, length(gamma_grid)),
    status = rep("not_required_after_global_infeasibility", length(gamma_grid)),
    feasible = rep(NA, length(gamma_grid)),
    solver = if (primal_reuse) "HiGHS-primal-simplex-dual-programme" else "HiGHS-dual-simplex", coordinate = k
  ))
  if (certificate$verified) {
    for (k in seq_len(p)) {
      grids[[k]]$status[] <- "feasible_by_same_full_matrix_certificate"
      grids[[k]]$feasible[] <- TRUE
    }
  }
  attempts <- list()
  cutoff <- NA_real_
  execution <- if (certificate$verified) which.min(gamma_grid) else order(gamma_grid, decreasing = TRUE)
  for (index in execution) {
    gamma <- gamma_grid[index]
    matrix_feasible <- TRUE
    for (k in seq_len(p)) {
      if (reset_basis) solver$clear_basis()
      solver$L(seq_len(2L * p), dantzig_dual_cost_v2(p, k, gamma))
      solver$solve(log_to_console = FALSE)
      status <- highs_model_status_label_v2(solver$status())
      grids[[k]]$status[index] <- status
      if (status %in% c("unbounded", "unbounded_or_infeasible")) {
        if (certificate$verified) {
          attempts[[length(attempts) + 1L]] <- data.frame(
            gamma = gamma, coordinate = k, status = status, residual = NA_real_, feasible = FALSE)
          condition <- simpleError("The iid CLIME dual status contradicts the verified primal feasibility certificate.")
          condition$failure_stage <- "precision_solver"
          condition$diagnostics <- list(gamma = gamma, coordinate = k, status = status,
            control = control, feasibility_certificate = certificate,
            solver_info = solver$info(), basis = solver$get_basis(), completed_cells = attempts, grids = grids)
          stop(condition)
        }
        grids[[k]]$feasible[index] <- FALSE
        matrix_feasible <- FALSE
        cutoff <- gamma
        attempts[[length(attempts) + 1L]] <- data.frame(
          gamma = gamma, coordinate = k, status = status, residual = NA_real_, feasible = FALSE
        )
        break
      }
      if (!identical(status, "optimal")) {
        attempts[[length(attempts) + 1L]] <- data.frame(
          gamma = gamma, coordinate = k, status = status, residual = NA_real_, feasible = FALSE)
        condition <- simpleError("The iid CLIME matrix solver returned an unaccepted numerical status.")
        condition$failure_stage <- "precision_solver"
        condition$diagnostics <- list(gamma = gamma, coordinate = k, status = status,
                                      control = control, feasibility_certificate = certificate,
                                      solver_info = solver$info(), basis = solver$get_basis(),
                                      completed_cells = attempts, grids = grids)
        stop(condition)
      }
      solution <- solver$solution()
      grids[[k]]$solver_records[[index]] <- list(
        info = solver$info(), basis = solver$get_basis(), solution = solution
      )
      omega <- dantzig_dual_omega_v2(solution$row_dual, p)
      canonical <- numeric(p)
      canonical[k] <- 1
      residual <- if (length(omega) == p && all(is.finite(omega)))
        max(abs(as.numeric(H %*% omega) - canonical)) else Inf
      grids[[k]]$omega[[index]] <- omega
      grids[[k]]$residual[index] <- residual
      grids[[k]]$l1_norm[index] <- sum(abs(omega))
      grids[[k]]$feasible[index] <- is.finite(residual) && residual <= gamma + 1e-6
      attempts[[length(attempts) + 1L]] <- data.frame(
        gamma = gamma, coordinate = k, status = status, residual = residual,
        feasible = is.finite(residual) && residual <= gamma + 1e-6
      )
      if (!is.finite(residual) || residual > gamma + 1e-6) {
        condition <- simpleError("The iid CLIME matrix solver failed its residual check.")
        condition$failure_stage <- "precision_solver"
        condition$diagnostics <- list(gamma = gamma, coordinate = k, status = status,
                                      residual = residual, omega = omega, control = control,
                                      feasibility_certificate = certificate,
                                      solver_info = solver$info(), basis = solver$get_basis(),
                                      completed_cells = attempts, grids = grids)
        stop(condition)
      }
    }
    if (!matrix_feasible) break
  }
  # 完整矩阵在某个gamma不可行时，所有更小gamma也不可行。
  list(grids = grids, attempts = if (length(attempts)) do.call(rbind, attempts) else data.frame(),
       global_infeasible_at = cutoff, model_constructions = 1L,
       solver_control = control,
       basis_policy = if (primal_reuse) "reuse_primal_feasible" else if (reset_basis) "clear_before_every_row" else "reuse_across_rows",
       feasibility_certificate = certificate,
       source_precision_definition = "common_gamma_full_matrix_CLIME")
}
