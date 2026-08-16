# Corrected reference implementation for the SJS reconstruction.
#
# This file is intentionally self-contained. It prioritises mathematical
# equivalence to docs/METHOD_SPECIFICATION.md over speed. Optimised code must
# reproduce this implementation on deterministic small and moderate examples.

`%||%` <- function(x, y) if (is.null(x)) y else x

assert_numeric_matrix_v2 <- function(x, name) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (any(!is.finite(x))) stop(name, " contains non-finite values.")
  x
}

validate_tau_h_v2 <- function(tau, h) {
  if (length(tau) != 1L || !is.finite(tau) || tau <= 0 || tau >= 1) {
    stop("tau must be a finite scalar in (0, 1).")
  }
  if (length(h) != 1L || !is.finite(h) || h <= 0) {
    stop("h must be a finite positive scalar.")
  }
  invisible(TRUE)
}

make_lambda_matrix_v2 <- function(lambda_gamma, q) {
  if (length(lambda_gamma) == 1L) {
    if (!is.finite(lambda_gamma) || lambda_gamma <= 0) {
      stop("lambda_gamma must be positive.")
    }
    Lambda <- diag(as.numeric(lambda_gamma), q)
  } else if (is.vector(lambda_gamma) && length(lambda_gamma) == q) {
    if (any(!is.finite(lambda_gamma)) || any(lambda_gamma <= 0)) {
      stop("All diagonal nuisance penalties must be positive.")
    }
    Lambda <- diag(as.numeric(lambda_gamma), q)
  } else {
    Lambda <- assert_numeric_matrix_v2(lambda_gamma, "lambda_gamma")
    if (!all(dim(Lambda) == c(q, q))) stop("lambda_gamma matrix has wrong dimension.")
    Lambda <- (Lambda + t(Lambda)) / 2
    ev <- eigen(Lambda, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) <= 0) stop("lambda_gamma matrix must be positive definite.")
  }
  Lambda
}

# -----------------------------------------------------------------------------
# Epanechnikov convolution smoothing
# -----------------------------------------------------------------------------

kernel_epan_pdf_v2 <- function(u) {
  out <- numeric(length(u))
  keep <- abs(u) <= 1
  out[keep] <- 0.75 * (1 - u[keep]^2)
  out
}

kernel_epan_cdf_v2 <- function(u) {
  out <- numeric(length(u))
  out[u >= 1] <- 1
  mid <- u > -1 & u < 1
  out[mid] <- 0.5 + 0.75 * (u[mid] - u[mid]^3 / 3)
  out
}

.epan_m1_primitive_v2 <- function(v) {
  (3 / 8) * v^2 - (3 / 16) * v^4
}

rho_smooth_v2 <- function(u, tau, h) {
  validate_tau_h_v2(tau, h)
  u <- as.numeric(u)
  x <- u / h
  cut <- pmin(1, pmax(-1, x))

  A1 <- kernel_epan_cdf_v2(cut)
  A2 <- 1 - A1
  m_minus <- .epan_m1_primitive_v2(-1)
  m_plus <- .epan_m1_primitive_v2(1)
  B1 <- .epan_m1_primitive_v2(cut) - m_minus
  B2 <- m_plus - .epan_m1_primitive_v2(cut)

  tau * (u * A1 - h * B1) + (tau - 1) * (u * A2 - h * B2)
}

psi_smooth_v2 <- function(u, tau, h) {
  validate_tau_h_v2(tau, h)
  tau - kernel_epan_cdf_v2(-as.numeric(u) / h)
}

phi_smooth_v2 <- function(u, tau, h) {
  validate_tau_h_v2(tau, h)
  kernel_epan_pdf_v2(as.numeric(u) / h) / h
}

check_loss_v2 <- function(u, tau) {
  u <- as.numeric(u)
  u * (tau - as.numeric(u < 0))
}

# -----------------------------------------------------------------------------
# Cluster handling and nuisance profiling
# -----------------------------------------------------------------------------

cluster_index_v2 <- function(cluster_id) {
  if (anyNA(cluster_id)) stop("cluster_id contains missing values.")
  f <- factor(cluster_id, levels = unique(cluster_id))
  rows <- split(seq_along(f), f, drop = TRUE)
  list(rows = rows, levels = names(rows), n = length(rows), factor = f)
}

solve_linear_pd_v2 <- function(A, B, name = "matrix") {
  A <- (A + t(A)) / 2
  ch <- tryCatch(chol(A), error = function(e) NULL)
  if (is.null(ch)) stop(name, " is not numerically positive definite.")
  backsolve(ch, forwardsolve(t(ch), B))
}

profile_cluster_gamma_v2 <- function(y_i,
                                     X_i,
                                     Z_i,
                                     beta,
                                     tau,
                                     h,
                                     Lambda,
                                     gamma_init = NULL,
                                     control = list()) {
  y_i <- as.numeric(y_i)
  X_i <- assert_numeric_matrix_v2(X_i, "X_i")
  Z_i <- assert_numeric_matrix_v2(Z_i, "Z_i")
  beta <- as.numeric(beta)
  m_i <- length(y_i)
  q <- ncol(Z_i)
  if (nrow(X_i) != m_i || nrow(Z_i) != m_i) stop("Cluster dimensions do not match.")
  if (ncol(X_i) != length(beta)) stop("beta has wrong length.")
  if (!all(dim(Lambda) == c(q, q))) stop("Lambda has wrong dimension.")

  gamma_init <- gamma_init %||% rep(0, q)
  gamma_init <- as.numeric(gamma_init)
  if (length(gamma_init) != q) stop("gamma_init has wrong length.")

  reltol <- control$reltol %||% 1e-11
  maxit <- control$maxit %||% 300L
  grad_tol <- control$grad_tol %||% 1e-8
  newton_maxit <- control$newton_maxit %||% 100L
  armijo <- control$armijo %||% 1e-4
  backtrack_factor <- control$backtrack_factor %||% 0.5
  max_backtrack <- control$max_backtrack %||% 50L
  min_step <- control$min_step %||% 1e-14

  objective <- function(gamma) {
    r <- y_i - as.numeric(X_i %*% beta) - as.numeric(Z_i %*% gamma)
    mean(rho_smooth_v2(r, tau, h)) +
      0.5 * drop(crossprod(gamma, Lambda %*% gamma))
  }
  gradient <- function(gamma) {
    r <- y_i - as.numeric(X_i %*% beta) - as.numeric(Z_i %*% gamma)
    as.numeric(-crossprod(Z_i, psi_smooth_v2(r, tau, h)) / m_i + Lambda %*% gamma)
  }
  hessian <- function(gamma) {
    r <- y_i - as.numeric(X_i %*% beta) - as.numeric(Z_i %*% gamma)
    w <- phi_smooth_v2(r, tau, h) / m_i
    crossprod(Z_i, Z_i * w) + Lambda
  }

  # BFGS supplies a stable starting point across the compact-support smoothing
  # regions. Because Lambda is positive definite, every nuisance Hessian is
  # positive definite; a damped Newton refinement then enforces the strict
  # stationarity tolerance needed by the profile-geometry audit.
  opt <- stats::optim(
    par = gamma_init,
    fn = objective,
    gr = gradient,
    method = "BFGS",
    control = list(reltol = reltol, maxit = maxit)
  )

  gamma <- as.numeric(opt$par)
  value <- objective(gamma)
  grad <- gradient(gamma)
  newton_iterations <- 0L
  newton_backtracks <- 0L
  newton_last_step <- NA_real_
  newton_line_search_failed <- FALSE

  if (all(is.finite(gamma)) && is.finite(value) && all(is.finite(grad))) {
    for (iteration in seq_len(newton_maxit)) {
      if (max(abs(grad)) <= grad_tol) break

      H <- hessian(gamma)
      direction <- tryCatch(
        -as.numeric(solve_linear_pd_v2(H, grad, name = "nuisance Hessian")),
        error = function(e) rep(NA_real_, q)
      )
      slope <- sum(grad * direction)
      if (any(!is.finite(direction)) || !is.finite(slope) || slope >= 0) {
        direction <- -grad
        slope <- -sum(grad^2)
      }

      step <- 1
      accepted <- FALSE
      for (bt in 0:max_backtrack) {
        candidate <- gamma + step * direction
        candidate_value <- objective(candidate)
        if (is.finite(candidate_value) &&
            candidate_value <= value + armijo * step * slope) {
          accepted <- TRUE
          break
        }
        step <- step * backtrack_factor
        if (step < min_step) break
      }

      newton_backtracks <- newton_backtracks + bt
      if (!accepted) {
        newton_line_search_failed <- TRUE
        break
      }

      gamma <- candidate
      value <- candidate_value
      grad <- gradient(gamma)
      newton_iterations <- iteration
      newton_last_step <- step
      if (any(!is.finite(grad))) break
    }
  }

  grad <- gradient(gamma)
  gradient_max <- max(abs(grad))
  converged <- all(is.finite(gamma)) && is.finite(value) &&
    all(is.finite(grad)) && gradient_max <= grad_tol

  list(
    gamma = gamma,
    objective = as.numeric(value),
    gradient = grad,
    gradient_max = gradient_max,
    converged = converged,
    optim_convergence = opt$convergence,
    optim_message = opt$message %||% "",
    function_evaluations = unname(opt$counts[["function"]] %||% NA_integer_),
    gradient_evaluations = unname(opt$counts[["gradient"]] %||% NA_integer_),
    newton_iterations = newton_iterations,
    newton_backtracks = newton_backtracks,
    newton_last_step = newton_last_step,
    newton_line_search_failed = newton_line_search_failed
  )
}

resolve_gamma_init_v2 <- function(gamma_init, cluster_name, index, q) {
  if (is.null(gamma_init)) return(rep(0, q))
  if (is.list(gamma_init)) {
    value <- gamma_init[[cluster_name]] %||% gamma_init[[index]]
    return(as.numeric(value %||% rep(0, q)))
  }
  if (is.matrix(gamma_init)) {
    if (ncol(gamma_init) != q) stop("gamma_init matrix has wrong q.")
    if (!is.null(rownames(gamma_init)) && cluster_name %in% rownames(gamma_init)) {
      return(as.numeric(gamma_init[cluster_name, ]))
    }
    if (nrow(gamma_init) >= index) return(as.numeric(gamma_init[index, ]))
  }
  stop("Unsupported gamma_init representation.")
}
