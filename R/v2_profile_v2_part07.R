# ---------------------------------------------------------------------------
# Adjudicated alternative solver for the frozen Dantzig LP: HiGHS.
#
# WHY THIS EXISTS: the mu-CV term solves the same frozen LP
#   min ||omega||_1  s.t.  ||H omega - e_k||_inf <= mu
# 4 folds x 8 mu = 32 times per coordinate, about 75-80% of replication time.
# CLARABEL is a conic interior-point code and the program is a pure linear
# program with 4p rows and 2p columns, so the natural question was whether a
# simplex code is cheaper. This file answered it.
#
# RESULT: the swap wins on well-conditioned constructions and loses on the real
# profile Hessian, which is what the pipeline solves. One row over the frozen
# grid, p=500, n=200, registered cell A02:
#
#   Hessian     condition   CLARABEL compact      HiGHS         HiGHS status
#   synthetic       6.04   4.9 s / 5.5 s      3.3 s / 2.6 s    optimal at once
#   real         1.6e+19   15.8 s / 15.0 s    181.3 s / 181.6 s  time limit x3
#
# The HiGHS runs carried a 60 s limit per solve and three of four candidates
# reached it. The real Hessian is numerically singular, smallest singular value
# 2.589e-19 against a largest of 4.153e+00. Both paths select the same mu, so
# the program, the grid and the decision rule are unaffected. See
# COMPUTATIONAL_OPTIMISATION_HIGHS_PATH_REPORT.md section 5.
#
# STATUS: instrumented evidence. Enabled only via
# JPDMEQR_DANTZIG_SOLVER=highs; the default remains the frozen CVXR path, the
# recommended server configuration remains clarabel_native with parallel mu-CV
# cores, and no run uses this switch.
#
# SCOPE: the LP handed to the solver has the same A, the same objective, the
# same row bounds, the same mu grid and the same smallest-feasible rule as in
# part06. Feasibility uses the primal residual against mu + 1e-6 and never the
# status label. Parity against the frozen CVXR path and the CLARABEL-native path
# is enforced by tests/testthat/test-dantzig-highs-parity.R.
#
# NOTE: the agreement with CLARABEL (7.4e-08 in omega on the constructions) is
# far below the acceptance threshold of the frozen rule, but it is not machine
# precision.
# ---------------------------------------------------------------------------

# HiGHS model-status code -> label. Order follows the HighsModelStatus enum:
# 0 NotSet, 1 LoadError, 2 ModelError, 3 PresolveError, 4 SolveError,
# 5 PostsolveError, 6 ModelEmpty, 7 Optimal, 8 Infeasible,
# 9 UnboundedOrInfeasible, 10 Unbounded, 11 ObjectiveBound, 12 ObjectiveTarget,
# 13 TimeLimit, 14 IterationLimit, 15 Unknown, 16 SolutionLimit, 17 Interrupt.
highs_model_status_label_v2 <- function(status) {
  s <- suppressWarnings(as.integer(status))
  if (length(s) != 1L || is.na(s)) return(paste0("status_", as.character(status)))
  labs <- c(
    "not_set", "load_error", "model_error", "presolve_error", "solve_error",
    "postsolve_error", "model_empty", "optimal", "infeasible",
    "unbounded_or_infeasible", "unbounded", "objective_bound",
    "objective_target", "time_limit", "iteration_limit", "unknown",
    "solution_limit", "interrupt"
  )
  if (s >= 0L && s < length(labs)) labs[s + 1L] else paste0("status_", s)
}

# Compact-form Dantzig row solved by HiGHS. Returns exactly the structure of
# solve_dantzig_row_v2 and solve_dantzig_row_native_v2, with the same
# ascending-mu / smallest-feasible semantics.
solve_dantzig_row_highs_v2 <- function(H,
                                       coordinate,
                                       mu_grid,
                                       solver_opts = list()) {
  if (!requireNamespace("highs", quietly = TRUE)) {
    stop("Package 'highs' is required for the HiGHS Dantzig path.")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required for the HiGHS Dantzig path.")
  }

  H <- assert_numeric_matrix_v2(H, "H")
  H <- (H + t(H)) / 2
  p <- nrow(H)
  if (ncol(H) != p) stop("H must be square.")

  if (is.character(coordinate)) {
    if (is.null(colnames(H)) || !coordinate %in% colnames(H)) {
      stop("Unknown coordinate name: ", coordinate)
    }
    k <- match(coordinate, colnames(H))
  } else {
    k <- as.integer(coordinate)
    if (length(k) != 1L || k < 1L || k > p) stop("Invalid coordinate index.")
  }

  mu_grid <- sort(unique(as.numeric(mu_grid)))
  if (length(mu_grid) == 0L || any(!is.finite(mu_grid)) || any(mu_grid <= 0)) {
    stop("mu_grid must contain positive finite values.")
  }

  e <- numeric(p)
  e[k] <- 1

  A <- dantzig_lp_compact_build_v2(H, p)
  L <- c(numeric(p), rep(1, p))
  lower <- c(rep(-Inf, p), rep(0, p))
  upper <- rep(Inf, 2L * p)
  lhs <- rep(-Inf, nrow(A))

  opts <- solver_opts
  if (is.null(opts$control)) opts$control <- highs::highs_control(threads = 1L)

  attempts <- vector("list", length(mu_grid))
  for (ii in seq_along(mu_grid)) {
    mu <- mu_grid[ii]
    b <- dantzig_lp_compact_rhs_v2(p, e, mu)
    res <- tryCatch(
      do.call(highs::highs_solve,
              c(list(L = L, lower = lower, upper = upper, A = A, lhs = lhs,
                     rhs = b), opts)),
      error = function(err) err
    )
    if (inherits(res, "error")) {
      attempts[[ii]] <- data.frame(mu = mu, status = "solver_error",
                                   residual = NA_real_,
                                   message = conditionMessage(res))
      next
    }
    status <- highs_model_status_label_v2(res$status)
    omega <- if (is.null(res$primal_solution)) NULL else
      as.numeric(res$primal_solution[seq_len(p)])
    if (is.null(omega) || length(omega) != p || any(!is.finite(omega))) {
      attempts[[ii]] <- data.frame(mu = mu, status = status, residual = NA_real_,
                                   message = res$status_message %||% "")
      next
    }
    # Same acceptance criterion as the frozen CVXR path and the compact native
    # path (primal residual against mu + 1e-6).
    residual <- max(abs(as.numeric(H %*% omega) - e))
    attempts[[ii]] <- data.frame(mu = mu, status = status, residual = residual,
                                 message = if (identical(status, "optimal")) "" else
                                   (res$status_message %||% ""))
    if (residual <= mu + 1e-6) {
      return(list(
        omega = omega,
        coordinate = k,
        coordinate_name = colnames(H)[k] %||% as.character(k),
        mu = mu,
        residual = residual,
        solver = "HiGHS",
        solver_status = status,
        l1_norm = sum(abs(omega)),
        l2_norm = sqrt(sum(omega^2)),
        attempts = do.call(rbind, attempts[seq_len(ii)]),
        feasible = TRUE
      ))
    }
  }

  list(
    omega = rep(NA_real_, p),
    coordinate = k,
    coordinate_name = colnames(H)[k] %||% as.character(k),
    mu = NA_real_,
    residual = NA_real_,
    solver = "HiGHS",
    solver_status = "infeasible_on_grid",
    l1_norm = NA_real_,
    l2_norm = NA_real_,
    attempts = do.call(rbind, attempts),
    feasible = FALSE
  )
}