submission_clime_qp_candidate <- function(H, coordinate, gamma) {
  p <- nrow(H)
  canonical <- numeric(p)
  canonical[coordinate] <- 1
  selected <- union(coordinate, head(order(abs(H[, coordinate]), decreasing = TRUE), min(16L, p)))
  history <- list()
  previous <- numeric(p)
  for (iteration in seq_len(p)) {
    count <- length(selected)
    identity <- Matrix::Diagonal(count)
    P <- Matrix::bdiag(Matrix::Matrix(H[selected, selected, drop = FALSE], sparse = TRUE), Matrix::Diagonal(count, 0))
    A <- rbind(cbind(identity, -identity), cbind(-identity, -identity))
    model <- osqp::osqp(P = P, q = c(-canonical[selected], rep(gamma, count)), A = A,
      l = rep(-Inf, 2L * count), u = rep(0, 2L * count),
      pars = osqp::osqpSettings(verbose = FALSE, max_iter = 10000L, eps_abs = 1e-7, eps_rel = 1e-7,
        eps_prim_inf = 1e-9, eps_dual_inf = 1e-9, adaptive_rho_interval = 25L))
    model$WarmStart(x = c(previous[selected], abs(previous[selected])))
    fitted <- model$Solve()
    status <- fitted$info$status
    history[[iteration]] <- list(selected = selected, info = fitted$info)
    if (identical(status, "dual infeasible")) {
      # osqp 0.6.3.3的prim_inf_cert字段保存native delta_x。
      native_ray <- if (packageVersion("osqp") < "1.0.0") fitted$prim_inf_cert else fitted$dual_inf_cert
      ray <- numeric(p)
      ray[selected] <- native_ray[seq_len(count)]
      ray <- ray / sum(abs(ray))
      residual <- max(abs(H %*% ray))
      margin <- ray[coordinate] - gamma
      if (any(!is.finite(ray)) || residual > 1e-7 * max(1, max(abs(H))) || margin <= 1e-7)
        stop("The QP ray failed the full CLIME infeasibility conditions.")
      return(list(status = "infeasible", qp_history = history, ray_residual = residual, ray_margin = margin))
    }
    if (!status %in% c("solved", "solved inaccurate")) {
      condition <- simpleError("The CLIME quadratic working problem did not converge.")
      condition$failure_stage <- "precision_solver"
      condition$diagnostics <- list(coordinate = coordinate, gamma = gamma, status = status, history = history)
      stop(condition)
    }
    omega <- numeric(p)
    omega[selected] <- fitted$x[seq_len(count)]
    if (any(!is.finite(omega))) stop("The CLIME quadratic candidate is nonfinite.")
    gradient <- as.numeric(H[, selected, drop = FALSE] %*% omega[selected] - canonical)
    omitted <- setdiff(seq_len(p), selected)
    violations <- omitted[abs(gradient[omitted]) > gamma + 1e-6]
    if (!length(violations)) return(list(status = "candidate", omega = omega, qp_history = history))
    selected <- union(selected, head(violations[order(abs(gradient[violations]), decreasing = TRUE)], 16L))
    previous <- omega
  }
  stop("The CLIME quadratic working set exhausted the coordinate budget.")
}
