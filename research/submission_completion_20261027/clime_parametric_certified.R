submission_clime_parametric_candidate <- function(H, coordinate, gamma) {
  invisible(loadNamespace("fastclime"))
  if (packageVersion("fastclime") != "1.4.1.1") stop("The certified parametric interface requires fastclime 1.4.1.1.")
  p <- nrow(H)
  canonical <- numeric(p)
  canonical[coordinate] <- 1
  maximum <- 4L * p
  limits <- unique(pmin(maximum, as.integer(128 * 2^(0:max(0, ceiling(log2(maximum / 128)))))))
  for (steps in limits) {
    path <- .C("dantzig", as.double(H), as.double(canonical), beta = numeric(p * steps),
      as.integer(p), as.double(gamma), as.integer(steps), gamma = rep(NA_real_, steps),
      NAOK = TRUE, PACKAGE = "fastclime")
    indices <- which(!is.na(path$gamma))
    knots <- path$gamma[indices]
    if (!length(knots)) stop("The native CLIME path did not return a knot.")
    if (any(is.finite(knots) & knots <= gamma)) {
      lower <- which(is.finite(knots) & knots <= gamma)[1L]
      upper <- max(1L, lower - 1L)
      coefficients <- matrix(path$beta, p, steps)
      if (knots[lower] == gamma || upper == lower) omega <- coefficients[, indices[lower]] else {
        weight <- (gamma - knots[lower]) / (knots[upper] - knots[lower])
        omega <- weight * coefficients[, indices[upper]] + (1 - weight) * coefficients[, indices[lower]]
      }
      if (any(!is.finite(omega))) stop("The interpolated CLIME path is nonfinite.")
      return(list(omega = omega, status = "candidate", knots = length(indices), budget = steps))
    }
    if (length(indices) < steps && all(is.finite(knots)) && tail(knots, 1L) > gamma)
      return(list(status = "infeasible", knots = length(indices), budget = steps))
  }
  condition <- simpleError("The native CLIME path exhausted its declared pivot budget.")
  condition$failure_stage <- "precision_solver"
  condition$diagnostics <- list(coordinate = coordinate, gamma = gamma, status = "parametric_iteration_limit",
    maximum_pivots = maximum, final_knot = tail(knots, 1L))
  stop(condition)
}

submission_clime_parametric_certified_grid <- function(H, gamma_grid, solver_opts = list()) {
  p <- nrow(H)
  gamma_grid <- sort(unique(gamma_grid))
  control <- highs::highs_control(threads = 1L, solver = "simplex", presolve = "on", log_to_console = FALSE,
    simplex_strategy = 1L, simplex_scale_strategy = 2L, simplex_iteration_limit = 20L * p,
    primal_feasibility_tolerance = 1e-8, dual_feasibility_tolerance = 1e-8)
  grids <- lapply(seq_len(p), function(k) list(mu = gamma_grid, omega = vector("list", length(gamma_grid)),
    solver_records = vector("list", length(gamma_grid)), residual = rep(NA_real_, length(gamma_grid)),
    l1_norm = rep(NA_real_, length(gamma_grid)), status = rep("not_required_after_global_infeasibility", length(gamma_grid)),
    feasible = rep(NA, length(gamma_grid)), solver = "CLIME_full_primal_dual_certificate", coordinate = k))
  records <- list()
  models <- 0L
  cutoff <- NA_real_
  for (index in order(gamma_grid, decreasing = TRUE)) {
    gamma <- gamma_grid[index]
    for (k in seq_len(p)) {
      canonical <- numeric(p)
      canonical[k] <- 1
      certificate <- submission_clime_certificate(H, k, gamma)
      native <- list(status = "one_coordinate_certificate")
      history <- list()
      status <- "optimal"
      if (!certificate$verified) {
        use_qp <- identical(solver_opts$initializer, "qp") ||
          (identical(solver_opts$initializer, "hybrid") && gamma >= .25)
        native <- if (use_qp) submission_clime_qp_candidate(H, k, gamma) else
          submission_clime_parametric_candidate(H, k, gamma)
        if (native$status == "infeasible") status <- "infeasible" else {
          certificate <- submission_clime_certificate(H, k, gamma, native$omega)
          if (!certificate$verified) {
            defect <- drop(H %*% native$omega) - canonical
            rows <- union(k, which(abs(abs(defect) - gamma) <= 1e-4))
            for (iteration in seq_len(p)) {
              block <- H[rows, , drop = FALSE]
              model <- highs::highs_model(L = rep(1, 2L * p), lower = rep(0, 2L * p), upper = rep(Inf, 2L * p),
                A = cbind(block, -block), lhs = canonical[rows] - gamma, rhs = canonical[rows] + gamma)
              solver <- highs::highs_solver(model, control)
              solver$solve(log_to_console = FALSE)
              status <- highs_model_status_label_v2(solver$status())
              history[[iteration]] <- list(rows = rows, status = status, info = solver$info())
              models <- models + 1L
              if (status == "infeasible") break
              if (status != "optimal") {
                condition <- simpleError("The CLIME constraint refinement did not reach an accepted status.")
                condition$failure_stage <- "precision_solver"
                condition$diagnostics <- list(gamma = gamma, coordinate = k, status = status, history = history)
                stop(condition)
              }
              solution <- solver$solution()
              omega <- solution$col_value[seq_len(p)] - solution$col_value[p + seq_len(p)]
              dual <- numeric(p)
              dual[rows] <- solution$row_dual
              defect <- drop(H %*% omega) - canonical
              objective <- sum(abs(omega))
              gap <- objective - dual[k] + gamma * sum(abs(dual))
              dual_residual <- max(abs(crossprod(H, dual)))
              certificate <- list(verified = max(abs(defect)) <= gamma + 1e-7 && dual_residual <= 1 + 1e-7 &&
                abs(gap) <= 1e-7 * max(1, objective), omega = omega, dual = dual,
                residual = max(abs(defect)), dual_residual = dual_residual, dual_gap = gap, objective = objective)
              if (certificate$verified) break
              violations <- setdiff(which(abs(defect) > gamma + 1e-7), rows)
              if (!length(violations)) stop("The CLIME refinement failed its full dual certificate.")
              rows <- union(rows, head(violations[order(abs(defect[violations]), decreasing = TRUE)], 16L))
            }
          }
        }
      }
      grids[[k]]$status[index] <- status
      if (status == "infeasible") {
        grids[[k]]$feasible[index] <- FALSE
        cutoff <- gamma
        records[[length(records) + 1L]] <- data.frame(gamma = gamma, coordinate = k, status = status,
          residual = NA_real_, feasible = FALSE, dual_gap = NA_real_, refinements = length(history))
        break
      }
      if (!certificate$verified || any(!is.finite(unlist(certificate)))) stop("The CLIME full matrix certificate was not attained.")
      grids[[k]]$omega[[index]] <- certificate$omega
      grids[[k]]$solver_records[[index]] <- list(native = native[c("status", "knots", "budget", "qp_history")], history = history,
        dual = certificate$dual, dual_gap = certificate$dual_gap, dual_residual = certificate$dual_residual)
      grids[[k]]$residual[index] <- certificate$residual
      grids[[k]]$l1_norm[index] <- certificate$objective
      grids[[k]]$feasible[index] <- TRUE
      records[[length(records) + 1L]] <- data.frame(gamma = gamma, coordinate = k, status = "optimal",
        residual = certificate$residual, feasible = TRUE, dual_gap = certificate$dual_gap, refinements = length(history))
    }
    if (is.finite(cutoff)) break
  }
  list(grids = grids, attempts = do.call(rbind, records), global_infeasible_at = cutoff,
    model_constructions = models, solver_control = control, basis_policy = "native_path_full_primal_dual_certificate",
    feasibility_certificate = list(verified = FALSE, definition = "per_row_primal_dual_certificates"),
    source_precision_definition = "common_gamma_full_matrix_CLIME")
}
