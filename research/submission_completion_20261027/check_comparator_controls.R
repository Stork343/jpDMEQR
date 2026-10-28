source("research/submission_completion_20261027/load_comparators.R")
records <- list()
for (tau in c(.25, .5, .75)) {
  set.seed(202610280L + as.integer(tau * 100))
  nobs <- sample(3:8, 40L, replace = TRUE)
  X <- cbind("(Intercept)" = 1, matrix(rnorm(sum(nobs) * 15L), sum(nobs)))
  colnames(X)[-1L] <- paste0("x", 1:15)
  attr(X, "assign") <- c(0L, 1:15)
  beta <- c(.3, 1, -.6, rep(0, 13L))
  y <- drop(X %*% beta) + rep(rnorm(40L), nobs) + rnorm(sum(nobs)) - qnorm(tau)
  control <- geeVerse::qpgeeControl(standardize = FALSE, maxit = 1000L, epsilon = 1e-4)
  reference <- geeVerse::qpgee(X, y, nobs, tau, control = control)
  candidate <- submission_qpgee(X, y, nobs, tau, control = control)
  difference <- max(abs(reference$coefficients - candidate$coefficients))
  if (difference > 1e-6 || reference$converged != candidate$converged ||
      abs(reference$hbic - candidate$hbic) > 1e-7) stop("The complete QGEE HBIC procedure changed.")
  records[[length(records) + 1L]] <- data.frame(method = "QGEE-SCAD", tau = tau,
    difference = difference, converged = candidate$converged)
  train <- list(y = y, X = X[, -1L], Z = matrix(1, length(y), 1L), cluster_id = rep(seq_len(40L), nobs))
  fitted <- fit_benchmark_double_pen_qlmm_v2(train, tau, c("x1", "x2"), list(), 20261028L,
    control = list(raw_sum_penalty_scales = TRUE, lambda_beta_grid = c(.25, .5, 1, 2), lambda_alpha_grid = c(.5, 1)))
  if (is.null(fitted$fit_object)) stop("The declared double-penalty programme did not return a fitted object.")
  expected <- c(beta = sqrt(tau * (1 - tau) * length(y) * log(15)),
    alpha = sqrt(tau * (1 - tau) * mean(nobs) * log(40)))
  scale_difference <- max(abs(fitted$fit_object$penalty_scales - expected))
  if (scale_difference > 1e-12) stop("The double-penalty scales do not match their specification.")
  source_fit <- fitted$fit_object$fit
  if (!all(source_fit$lambda_beta_grid %in% (c(.25, .5, 1, 2) * expected["beta"])) ||
      !all(source_fit$lambda_alpha_grid %in% (c(.5, 1) * expected["alpha"]))) stop("The double-penalty grid is incomplete.")
  records[[length(records) + 1L]] <- data.frame(method = "DOUBLE-PEN-QLMM", tau = tau,
    difference = scale_difference, converged = fitted$converged)
}
write.csv(do.call(rbind, records), "results/scratch/submission_completion_20261027/comparator_controls_checks.csv", row.names = FALSE)
print(do.call(rbind, records))
