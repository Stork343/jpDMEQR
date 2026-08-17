fit_profile_lasso_v2 <- function(y,
                                 X,
                                 Z,
                                 cluster_id,
                                 tau,
                                 h,
                                 lambda_beta,
                                 lambda_gamma = 1,
                                 penalty_factor = NULL,
                                 beta_init = NULL,
                                 control = list()) {
  validate_tau_h_v2(tau, h)
  y <- as.numeric(y)
  X <- assert_numeric_matrix_v2(X, "X")
  Z <- assert_numeric_matrix_v2(Z, "Z")
  p <- ncol(X)
  if (!is.finite(lambda_beta) || lambda_beta < 0) stop("lambda_beta must be nonnegative.")

  penalty_factor <- penalty_factor %||% rep(1, p)
  penalty_factor <- as.numeric(penalty_factor)
  if (length(penalty_factor) != p || any(!is.finite(penalty_factor)) ||
      any(penalty_factor < 0)) stop("Invalid penalty_factor.")

  beta <- as.numeric(beta_init %||% rep(0, p))
  if (length(beta) != p) stop("beta_init has wrong length.")

  max_iter <- control$max_iter %||% 200L
  beta_tol <- control$beta_tol %||% 1e-6
  kkt_tol <- control$kkt_tol %||% 1e-5
  backtrack_factor <- control$backtrack_factor %||% 0.5
  max_backtrack <- control$max_backtrack %||% 30L
  min_step <- control$min_step %||% 1e-12
  verbose <- isTRUE(control$verbose)
  nuisance_control <- control$nuisance_control %||% list()

  current <- profile_components_v2(
    y, X, Z, cluster_id, beta, tau, h, lambda_gamma,
    gamma_init = NULL, need_hessian = TRUE,
    nuisance_control = nuisance_control
  )
  composite <- current$objective + lambda_beta * sum(penalty_factor * abs(beta))
  history <- vector("list", max_iter)
  converged <- FALSE
  failure_stage <- ""

  for (iter in seq_len(max_iter)) {
    eigmax <- largest_eigenvalue_power_v2(current$hessian)
    step <- 1 / max(eigmax, 1e-6)
    score <- current$score
    accepted <- FALSE
    candidate <- NULL

    for (bt in 0:max_backtrack) {
      beta_new <- soft_threshold_weighted_v2(
        beta - step * score,
        step * lambda_beta * penalty_factor
      )
      delta <- beta_new - beta
      candidate <- tryCatch(
        profile_components_v2(
          y, X, Z, cluster_id, beta_new, tau, h, lambda_gamma,
          gamma_init = current$gamma, need_hessian = TRUE,
          nuisance_control = nuisance_control
        ),
        error = function(e) e
      )
      if (!inherits(candidate, "error")) {
        upper <- current$objective + sum(score * delta) +
          sum(delta^2) / (2 * step) + 1e-12
        if (candidate$objective <= upper) {
          accepted <- TRUE
          break
        }
      }
      step <- step * backtrack_factor
      if (step < min_step) break
    }

    if (!accepted) {
      failure_stage <- "profile_backtracking"
      break
    }

    composite_new <- candidate$objective +
      lambda_beta * sum(penalty_factor * abs(beta_new))
    beta_change <- max(abs(beta_new - beta))
    kkt <- kkt_residual_profile_v2(
      beta_new, candidate$score, lambda_beta, penalty_factor
    )

    history[[iter]] <- data.frame(
      iteration = iter,
      profile_objective = candidate$objective,
      composite_objective = composite_new,
      beta_change = beta_change,
      kkt_residual = kkt,
      step = step,
      backtracks = bt,
      nuisance_converged = candidate$nuisance_converged,
      max_nuisance_gradient = candidate$max_nuisance_gradient,
      profile_identity_error = candidate$profile_identity_error,
      stringsAsFactors = FALSE
    )

    if (verbose) {
      message(sprintf(
        "iter=%d composite=%.8f change=%.3e kkt=%.3e step=%.3e",
        iter, composite_new, beta_change, kkt, step
      ))
    }

    beta <- beta_new
    current <- candidate
    composite <- composite_new

    if (beta_change <= beta_tol && kkt <= kkt_tol && current$nuisance_converged) {
      converged <- TRUE
      break
    }
  }

  history <- do.call(rbind, history[!vapply(history, is.null, logical(1))])
  final_kkt <- kkt_residual_profile_v2(
    beta, current$score, lambda_beta, penalty_factor
  )

  list(
    beta = setNames(beta, colnames(X) %||% paste0("x", seq_len(p))),
    gamma = current$gamma,
    components = current,
    profile_objective = current$objective,
    composite_objective = composite,
    lambda_beta = lambda_beta,
    lambda_gamma = lambda_gamma,
    penalty_factor = penalty_factor,
    converged = converged,
    failure_stage = failure_stage,
    kkt_residual = final_kkt,
    iterations = nrow(history),
    history = history,
    call = match.call()
  )
}

# -----------------------------------------------------------------------------
# Row-wise Dantzig direction and debiased inference
# -----------------------------------------------------------------------------

solve_dantzig_row_v2 <- function(H,
                                 coordinate,
                                 mu_grid,
                                 solver_preference = c("CLARABEL", "ECOS", "SCS"),
                                 solver_opts = list()) {
  if (!requireNamespace("CVXR", quietly = TRUE)) {
    stop("Package 'CVXR' is required for the Dantzig program.")
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

  installed <- CVXR::installed_solvers()
  solver <- solver_preference[solver_preference %in% installed][1]
  if (is.na(solver)) stop("No requested CVXR solver is installed.")

  e <- numeric(p)
  e[k] <- 1
  attempts <- vector("list", length(mu_grid))

  for (ii in seq_along(mu_grid)) {
    mu <- mu_grid[ii]
    omega_var <- CVXR::Variable(p)
    constraints <- list(
      H %*% omega_var - e <= mu,
      H %*% omega_var - e >= -mu
    )
    problem <- CVXR::Problem(
      CVXR::Minimize(CVXR::p_norm(omega_var, 1)),
      constraints
    )
    # CVXR >= 1.9 renamed the exported solver entry point from solve to
    # psolve; keep both spellings working across installed versions.
    cvxr_solve <- if ("psolve" %in% getNamespaceExports("CVXR")) {
      getExportedValue("CVXR", "psolve")
    } else {
      getExportedValue("CVXR", "solve")
    }
    solved <- tryCatch(
      do.call(cvxr_solve, c(list(problem = problem, solver = solver), solver_opts)),
      error = function(e) e
    )

    if (inherits(solved, "error")) {
      attempts[[ii]] <- data.frame(mu = mu, status = "solver_error",
                                   residual = NA_real_, message = conditionMessage(solved))
      next
    }

    status <- solved$status %||% "unknown"
    omega <- tryCatch(as.numeric(solved$getValue(omega_var)), error = function(e) NULL)
    if (is.null(omega) || length(omega) != p || any(!is.finite(omega))) {
      attempts[[ii]] <- data.frame(mu = mu, status = status,
                                   residual = NA_real_, message = "No finite solution returned")
      next
    }
    residual <- max(abs(as.numeric(H %*% omega - e)))
    attempts[[ii]] <- data.frame(mu = mu, status = status,
                                 residual = residual, message = "")

    if (residual <= mu + 1e-6 && status %in% c("optimal", "optimal_inaccurate")) {
      return(list(
        omega = omega,
        coordinate = k,
        coordinate_name = colnames(H)[k] %||% as.character(k),
        mu = mu,
        residual = residual,
        solver = solver,
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
    solver = solver,
    solver_status = "infeasible_on_grid",
    l1_norm = NA_real_,
    l2_norm = NA_real_,
    attempts = do.call(rbind, attempts),
    feasible = FALSE
  )
}
