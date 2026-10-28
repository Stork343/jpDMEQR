submission_clime_certificate <- function(H, coordinate, gamma, initial = NULL) {
  p <- nrow(H)
  canonical <- numeric(p)
  canonical[coordinate] <- 1
  omega <- dual <- numeric(p)
  if (is.null(initial)) {
    j <- which.max(abs(H[coordinate, ]))
    if (abs(H[coordinate, j]) < 1e-14) return(list(verified = FALSE))
    omega[j] <- (1 - gamma) / H[coordinate, j]
    dual[coordinate] <- 1 / abs(H[coordinate, j])
  } else {
    support <- which(abs(initial) > 1e-6)
    if (!length(support)) return(list(verified = FALSE))
    defect <- as.numeric(H[, support, drop = FALSE] %*% initial[support] - canonical)
    ordered <- order(abs(abs(defect) - gamma))
    decomposition <- qr(t(H[ordered, support, drop = FALSE]), tol = 1e-10, LAPACK = FALSE)
    if (decomposition$rank != length(support)) return(list(verified = FALSE))
    active <- ordered[decomposition$pivot[seq_along(support)]]
    block <- H[active, support, drop = FALSE]
    if (rcond(block) < 1e-10) return(list(verified = FALSE))
    omega[support] <- solve(block, canonical[active] + gamma * sign(defect[active]))
    dual[active] <- solve(t(block), sign(omega[support]))
  }
  support <- which(omega != 0)
  active <- which(dual != 0)
  residual <- max(abs(as.numeric(H[, support, drop = FALSE] %*% omega[support]) - canonical))
  dual_residual <- max(abs(as.numeric(crossprod(H[active, , drop = FALSE], dual[active]))))
  objective <- sum(abs(omega))
  gap <- objective - dual[coordinate] + gamma * sum(abs(dual))
  verified <- all(is.finite(c(omega, dual, residual, gap))) && residual <= gamma + 1e-7 &&
    dual_residual <= 1 + 1e-7 && abs(gap) <= 1e-7 * max(1, objective)
  list(verified = verified, omega = omega, dual = dual, residual = residual,
    dual_residual = dual_residual, dual_gap = gap, objective = objective)
}

submission_clime_certified_grid <- function(H, gamma_grid, solver_opts = list()) {
  p <- nrow(H)
  gamma_grid <- sort(unique(gamma_grid))
  control <- highs::highs_control(threads = 1L, solver = "simplex", presolve = "off", log_to_console = FALSE,
    simplex_strategy = 1L, simplex_scale_strategy = 2L, simplex_iteration_limit = 20L * p,
    primal_feasibility_tolerance = 1e-8, dual_feasibility_tolerance = 1e-8)
  model <- highs::highs_model(L = rep(1, 2L * p), lower = rep(0, 2L * p), upper = rep(Inf, 2L * p),
    A = cbind(H, -H), lhs = rep(-1, p), rhs = rep(1, p))
  solver <- highs::highs_solver(model, control)
  grids <- lapply(seq_len(p), function(k) list(mu = gamma_grid, omega = vector("list", length(gamma_grid)),
    solver_records = vector("list", length(gamma_grid)), residual = rep(NA_real_, length(gamma_grid)),
    l1_norm = rep(NA_real_, length(gamma_grid)), status = rep("not_required_after_global_infeasibility", length(gamma_grid)),
    feasible = rep(NA, length(gamma_grid)), solver = "certified_CLIME", coordinate = k))
  records <- list()
  cutoff <- NA_real_
  for (index in order(gamma_grid, decreasing = TRUE)) {
    gamma <- gamma_grid[index]
    warm_directions <- NULL
    for (k in seq_len(p)) {
      canonical <- numeric(p)
      canonical[k] <- 1
      certified <- submission_clime_certificate(H, k, gamma)
      origin <- "one_coordinate_certificate"
      if (!certified$verified && k > 1L) {
        if (is.null(warm_directions)) warm_directions <- submission_flare_matrix(H, gamma,
          precision = 1e-5, max_iterations = 1000L)
        initial <- warm_directions$matrices[[1L]][, k]
        certified <- submission_clime_certificate(H, k, gamma, initial)
        origin <- "active_constraints_certificate"
      }
      status <- "optimal"
      info <- list(simplex_iteration_count = 0L, certificate = origin)
      if (!certified$verified) {
        solver$cbounds(seq_len(p), canonical - gamma, canonical + gamma)
        solver$solve(log_to_console = FALSE)
        status <- highs_model_status_label_v2(solver$status())
        info <- solver$info()
        if (status == "infeasible") {
          grids[[k]]$status[index] <- status
          grids[[k]]$feasible[index] <- FALSE
          cutoff <- gamma
          records[[length(records) + 1L]] <- data.frame(gamma = gamma, coordinate = k, status = status,
            residual = NA_real_, feasible = FALSE, dual_gap = NA_real_, certificate = "HiGHS_infeasibility")
          break
        }
        if (status != "optimal") {
          condition <- simpleError("The full CLIME programme did not reach an accepted numerical status.")
          condition$failure_stage <- "precision_solver"
          condition$diagnostics <- list(gamma = gamma, coordinate = k, status = status,
            info = info, control = control, completed_cells = records)
          stop(condition)
        }
        solution <- solver$solution()
        omega <- solution$col_value[seq_len(p)] - solution$col_value[p + seq_len(p)]
        dual <- solution$row_dual
        certified <- list(omega = omega, dual = dual, residual = max(abs(H %*% omega - canonical)),
          dual_residual = max(abs(crossprod(H, dual))), objective = sum(abs(omega)),
          dual_gap = sum(abs(omega)) - dual[k] + gamma * sum(abs(dual)))
        if (any(!is.finite(unlist(certified))) || certified$residual > gamma + 1e-6 ||
            certified$dual_residual > 1 + 1e-6 || abs(certified$dual_gap) > 1e-6 * max(1, certified$objective))
          stop("The full CLIME programme failed its primal-dual certificate.")
        origin <- "HiGHS_primal_dual_certificate"
      }
      grids[[k]]$omega[[index]] <- certified$omega
      grids[[k]]$solver_records[[index]] <- list(info = info, certificate = origin, dual = certified$dual,
        dual_gap = certified$dual_gap, dual_residual = certified$dual_residual)
      grids[[k]]$residual[index] <- certified$residual
      grids[[k]]$l1_norm[index] <- certified$objective
      grids[[k]]$status[index] <- "optimal"
      grids[[k]]$feasible[index] <- TRUE
      records[[length(records) + 1L]] <- data.frame(gamma = gamma, coordinate = k, status = "optimal",
        residual = certified$residual, feasible = TRUE, dual_gap = certified$dual_gap, certificate = origin)
    }
    if (is.finite(cutoff)) break
  }
  list(grids = grids, attempts = do.call(rbind, records), global_infeasible_at = cutoff,
    model_constructions = 1L, solver_control = control, basis_policy = "certified_active_constraints",
    feasibility_certificate = list(verified = FALSE, definition = "per_row_primal_dual_certificates"),
    source_precision_definition = "common_gamma_full_matrix_CLIME")
}
