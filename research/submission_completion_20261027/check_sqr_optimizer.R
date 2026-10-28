source("research/submission_completion_20261027/load.R")
records <- list()
for (p in c(8L, 30L, 80L)) for (tau in c(.25, .5, .75)) {
  dat <- generate_profile_qr_data_v2(n_clusters = 50L, p = p, s = 5L, tau = tau,
    seed = 202610270L + p + as.integer(100 * tau))
  common <- list(max_iter = 5000L, beta_tol = 1e-10, kkt_tol = 1e-9)
  reference <- fit_sqr_iid_lasso_v2(dat$y, dat$X, tau, .4, .07, control = common)
  accelerated <- fit_sqr_iid_lasso_v2(dat$y, dat$X, tau, .4, .07,
    control = c(common, list(optimizer = "matrix_free_fista")))
  error <- max(abs(reference$beta - accelerated$beta))
  objective_error <- abs(reference$objective - accelerated$objective)
  if (!reference$converged || !accelerated$converged || error > 2e-5 || objective_error > 1e-9) stop("The Gaussian quantile solvers differ beyond numerical accuracy.")
  records[[length(records) + 1L]] <- data.frame(p = p, tau = tau, coefficient_error = error,
    objective_error = objective_error, reference_kkt = reference$kkt_residual,
    accelerated_kkt = accelerated$kkt_residual, reference_iterations = reference$iterations,
    accelerated_iterations = accelerated$iterations)
}
write.csv(do.call(rbind, records), "results/scratch/submission_completion_20261027/sqr_optimizer_checks.csv", row.names = FALSE)
print(do.call(rbind, records))
