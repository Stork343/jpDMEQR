source("research/submission_completion_20261027/load.R")
source("research/submission_completion_20261027/qgee_algebra.R")
Rcpp::sourceCpp("research/submission_completion_20261027/qgee_algebra.cpp",
  cacheDir = "results/scratch/submission_completion_20261027/qgee_cache", rebuild = FALSE)
namespace <- asNamespace("geeVerse")
reference_update <- get(".update_beta", namespace)
reference_kernel <- get(".fit_single_lambda", namespace)
candidate_kernel <- submission_qgee_kernel()
records <- list()
for (p in c(8L, 30L, 100L)) for (tau in c(.25, .5, .75)) {
  set.seed(202610280L + p + as.integer(tau * 100))
  nobs <- sample(1:8, 50L, replace = TRUE)
  X <- cbind(1, matrix(rnorm(sum(nobs) * (p - 1L)), sum(nobs)))
  colnames(X) <- c("(Intercept)", paste0("x", seq_len(p - 1L)))
  beta <- c(.2, 1, -.5, rep(0, p - 3L))
  y <- drop(X %*% beta) + rep(rnorm(50), nobs) + rnorm(sum(nobs))
  R <- matrix(.3, max(nobs), max(nobs))
  diag(R) <- 1
  arguments <- list(X, y, X %*% beta, beta, R, nobs, tau, 50L, c(0L, cumsum(nobs)), .1, c(0, rep(1, p - 1L)))
  reference <- do.call(reference_update, arguments)$beta_update
  candidate <- do.call(submission_qgee_update, arguments)$beta_update
  update_difference <- max(abs(reference - candidate))
  if (update_difference > 1e-8) stop("QGEE algebra changed its coefficient update.")
  control <- geeVerse::qpgeeControl(epsilon = 1e-4, maxit = 1000L, standardize = FALSE)
  initial <- as.numeric(quantreg::rq.fit.br(X, y, tau = tau)$coefficients)
  started <- proc.time()[[3L]]
  reference <- reference_kernel(X, y, nobs, tau, "exchangeable", .1, rep(1, p - 1L), TRUE, initial, control)
  reference_seconds <- proc.time()[[3L]] - started
  started <- proc.time()[[3L]]
  candidate <- candidate_kernel(X, y, nobs, tau, "exchangeable", .1, rep(1, p - 1L), TRUE, initial, control)
  candidate_seconds <- proc.time()[[3L]] - started
  coefficient_difference <- max(abs(reference$coefficients - candidate$coefficients))
  if (reference$converged != candidate$converged || coefficient_difference > 1e-6 ||
      abs(reference$hbic - candidate$hbic) > 1e-7 || !identical(reference$X_selected, candidate$X_selected))
    stop("QGEE fitted coefficients or selection changed after algebraic optimisation.")
  records[[length(records) + 1L]] <- data.frame(p = p, tau = tau, update_difference = update_difference,
    coefficient_difference = coefficient_difference, reference_seconds = reference_seconds, candidate_seconds = candidate_seconds,
    converged = candidate$converged, iterations = candidate$iterations)
  write.csv(do.call(rbind, records), "results/scratch/submission_completion_20261027/qgee_algebra_checks.csv", row.names = FALSE)
}
print(do.call(rbind, records))
