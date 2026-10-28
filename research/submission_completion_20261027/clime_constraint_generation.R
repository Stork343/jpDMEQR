submission_clime_constraint_grid <- function(H, gamma_grid, solver_opts = list()) {
  p <- nrow(H)
  gamma_grid <- sort(unique(gamma_grid))
  control <- if (identical(solver_opts$subproblem, "ipm")) highs::highs_control(
    threads = 1L, solver = "ipm", presolve = "on", log_to_console = FALSE, run_crossover = "off",
    ipm_iteration_limit = 200L, ipm_optimality_tolerance = 1e-10,
    primal_feasibility_tolerance = 1e-9, dual_feasibility_tolerance = 1e-9) else
    highs::highs_control(threads = 1L, solver = "simplex", presolve = "on", log_to_console = FALSE,
    simplex_strategy = 1L, simplex_scale_strategy = 2L, simplex_iteration_limit = 20L * p,
    primal_feasibility_tolerance = 1e-8, dual_feasibility_tolerance = 1e-8)
  grids <- lapply(seq_len(p), function(k) list(mu = gamma_grid, omega = vector("list", length(gamma_grid)),
    solver_records = vector("list", length(gamma_grid)), residual = rep(NA_real_, length(gamma_grid)),
    l1_norm = rep(NA_real_, length(gamma_grid)), status = rep("not_required_after_global_infeasibility", length(gamma_grid)),
    feasible = rep(NA, length(gamma_grid)), solver = "HiGHS_certified_constraint_generation", coordinate = k))
  records <- list()
  cutoff <- NA_real_
  models <- 0L
  for (index in order(gamma_grid, decreasing = TRUE)) {
    gamma <- gamma_grid[index]
    for (k in seq_len(p)) {
      canonical <- numeric(p)
      canonical[k] <- 1
      certified <- submission_clime_certificate(H, k, gamma)
      rows <- columns <- k
      history <- list()
      status <- "optimal"
      if (!certified$verified) for (iteration in seq_len(2L * p)) {
        block <- H[rows, columns, drop = FALSE]
        count <- length(columns)
        if (identical(solver_opts$subproblem, "clarabel")) {
          inequalities <- rbind(cbind(block, -block), cbind(-block, block), -diag(2L * count))
          fitted <- clarabel::clarabel(A = Matrix::Matrix(inequalities, sparse = TRUE),
            b = c(canonical[rows] + gamma, -canonical[rows] + gamma, rep(0, 2L * count)),
            q = rep(1, 2L * count), cones = list(l = nrow(inequalities)),
            control = list(verbose = FALSE, max_iter = 200L, tol_gap_abs = 1e-9,
              tol_gap_rel = 1e-9, tol_feas = 1e-9))
          status <- if (fitted$status %in% c(2L, 5L)) "optimal" else if (fitted$status %in% c(3L, 6L)) "infeasible" else
            paste0("clarabel_status_", fitted$status)
          solution <- list(col_value = fitted$x,
            row_dual = fitted$z[length(rows) + seq_along(rows)] - fitted$z[seq_along(rows)])
          info <- list(status = fitted$status, iterations = fitted$iterations,
            primal_objective = fitted$obj_val, dual_objective = fitted$obj_val_dual)
          if (status == "infeasible") {
            ray <- numeric(p)
            ray[rows] <- solution$row_dual
            ray <- ray / sum(abs(ray))
            ray_residual <- max(abs(crossprod(H, ray)))
            ray_margin <- ray[k] - gamma
            if (any(!is.finite(ray)) || ray_residual > 1e-8 * max(1, max(abs(H))) || ray_margin <= 1e-7)
              stop("The CLIME infeasibility certificate did not satisfy the full matrix conditions.")
            info$infeasibility_residual <- ray_residual
            info$infeasibility_margin <- ray_margin
          }
        } else {
          model <- highs::highs_model(L = rep(1, 2L * count), lower = rep(0, 2L * count),
            upper = rep(Inf, 2L * count), A = cbind(block, -block),
            lhs = canonical[rows] - gamma, rhs = canonical[rows] + gamma)
          solver <- highs::highs_solver(model, control)
          solver$solve(log_to_console = FALSE)
          status <- highs_model_status_label_v2(solver$status())
          info <- solver$info()
          solution <- solver$solution()
        }
        models <- models + 1L
        history[[iteration]] <- list(rows = rows, columns = columns, status = status, info = info)
        # 对PSD矩阵，H[rows, rows]和H[rows, ]具有相同的列空间。
        if (status == "infeasible") break
        if (status != "optimal") {
          condition <- simpleError("The restricted CLIME programme did not reach an accepted numerical status.")
          condition$failure_stage <- "precision_solver"
          condition$diagnostics <- list(gamma = gamma, coordinate = k, status = status, history = history, control = control)
          stop(condition)
        }
        omega <- dual <- numeric(p)
        omega[columns] <- solution$col_value[seq_len(count)] - solution$col_value[count + seq_len(count)]
        dual[rows] <- solution$row_dual
        defect <- as.numeric(H[, columns, drop = FALSE] %*% omega[columns] - canonical)
        pricing <- as.numeric(crossprod(H[rows, , drop = FALSE], dual[rows]))
        objective <- sum(abs(omega))
        gap <- objective - dual[k] + gamma * sum(abs(dual))
        certified <- list(verified = max(abs(defect)) <= gamma + 1e-7 && max(abs(pricing)) <= 1 + 1e-7 &&
          abs(gap) <= 1e-7 * max(1, objective), omega = omega, dual = dual,
          residual = max(abs(defect)), dual_residual = max(abs(pricing)), dual_gap = gap, objective = objective)
        if (certified$verified) break
        add_rows <- setdiff(which(abs(defect) > gamma + 1e-7), rows)
        add_columns <- setdiff(which(abs(pricing) > 1 + 1e-7), columns)
        if (!length(add_rows) && !length(add_columns)) stop("The CLIME restricted optimum failed its gap certificate.")
        add_rows <- head(add_rows[order(abs(defect[add_rows]), decreasing = TRUE)], 16L)
        add_columns <- head(add_columns[order(abs(pricing[add_columns]), decreasing = TRUE)], 16L)
        rows <- union(rows, add_rows)
        columns <- union(columns, c(add_rows, add_columns))
      }
      grids[[k]]$status[index] <- status
      if (status == "infeasible") {
        grids[[k]]$feasible[index] <- FALSE
        cutoff <- gamma
        records[[length(records) + 1L]] <- data.frame(gamma = gamma, coordinate = k, status = status,
          residual = NA_real_, feasible = FALSE, dual_gap = NA_real_, programmes = length(history))
        break
      }
      if (!certified$verified || any(!is.finite(unlist(certified)))) stop("The full CLIME certificate was not attained.")
      grids[[k]]$omega[[index]] <- certified$omega
      grids[[k]]$solver_records[[index]] <- list(history = history, dual = certified$dual,
        dual_gap = certified$dual_gap, dual_residual = certified$dual_residual)
      grids[[k]]$residual[index] <- certified$residual
      grids[[k]]$l1_norm[index] <- certified$objective
      grids[[k]]$feasible[index] <- TRUE
      records[[length(records) + 1L]] <- data.frame(gamma = gamma, coordinate = k, status = "optimal",
        residual = certified$residual, feasible = TRUE, dual_gap = certified$dual_gap, programmes = length(history))
    }
    if (is.finite(cutoff)) break
  }
  list(grids = grids, attempts = do.call(rbind, records), global_infeasible_at = cutoff,
    model_constructions = models, solver_control = control, basis_policy = "full_primal_and_dual_pricing",
    feasibility_certificate = list(verified = FALSE, definition = "per_row_primal_dual_certificates"),
    source_precision_definition = "common_gamma_full_matrix_CLIME")
}
