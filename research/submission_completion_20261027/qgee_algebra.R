submission_qgee_update <- function(x_active, y, mu, beta_current, R, nobs, tau, n_subjects,
                                  cluster_indices, lambda, penalty_factor_active, f0 = NULL,
                                  inverse_apply = NULL) {
  if (is.null(f0)) f0 <- rep(1, length(y))
  if (any(!is.finite(R))) stop("The QGEE working correlation is nonfinite.")
  namespace <- asNamespace("geeVerse")
  inverse_function <- get("geninv", namespace)
  penalty_derivative <- get("pp_scad_lin", namespace)
  radii <- sqrt(rowSums(x_active^2))
  radii[radii < .Machine$double.eps] <- 1e-6
  residual <- as.numeric(y - mu)
  argument <- sqrt(n_subjects) * residual / radii
  score <- pnorm(argument) - (1 - tau)
  density <- dnorm(argument) / radii
  inverses <- vector("list", max(nobs))
  for (m in unique(nobs[nobs > 0])) inverses[[m]] <- inverse_function(R[seq_len(m), seq_len(m), drop = FALSE])
  geometry <- submission_qgee_geometry_cpp(x_active, score, density, as.numeric(f0), nobs, inverses)
  hessian <- sqrt(n_subjects) * geometry$hessian
  score_vector <- geometry$score
  if (lambda > 0) {
    weights <- penalty_derivative(abs(as.vector(beta_current)), lambda) /
      (abs(as.vector(beta_current)) + 1e-6) * penalty_factor_active
    sigma <- n_subjects * diag(weights, nrow = ncol(x_active))
    hessian <- hessian + sigma
    score_vector <- score_vector - sigma %*% beta_current
  }
  update <- if (is.null(inverse_apply)) inverse_function(hessian) %*% score_vector else
    inverse_apply(hessian, matrix(score_vector, ncol = 1L))
  list(beta_update = update)
}

submission_qgee_kernel <- function(inverse_apply = submission_qgee_inverse_lapack) {
  namespace <- asNamespace("geeVerse")
  kernel <- get(".fit_single_lambda", namespace)
  accelerated <- new.env(parent = namespace)
  accelerated$.update_beta <- function(...) submission_qgee_update(..., inverse_apply = inverse_apply)
  environment(kernel) <- accelerated
  kernel
}

submission_qpgee <- function(x, y, nobs, tau = .5, ...) {
  arguments <- c(list(x = x, y = y, nobs = nobs, tau = tau), list(...))
  if (is.null(arguments$betaint) && isFALSE(arguments$control$standardize))
    arguments$betaint <- as.numeric(quantreg::rq.fit.br(x, y, tau = tau)$coefficients)
  method <- get("qpgee.default", asNamespace("geeVerse"))
  accelerated <- new.env(parent = environment(method))
  accelerated$.fit_single_lambda <- submission_qgee_kernel()
  environment(method) <- accelerated
  do.call(method, arguments)
}
