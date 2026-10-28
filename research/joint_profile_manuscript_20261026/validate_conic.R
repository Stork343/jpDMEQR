source("research/score_validity_20261026/load_joint.R")
source("research/joint_profile_manuscript_20261026/conic_search.R")
records <- list()
for (q in 1:2) {
  for (tau in c(.25, .5, .75)) {
    for (lambda in c(.5, 1, 2)) {
      dat <- generate_profile_qr_data_v2(n_clusters = 40L, p = 6L, s = 4L, tau = tau,
                                         q = q, seed = 202610270L + 100L * q + 10L * as.integer(4 * tau) + as.integer(2 * lambda))
      if (q == 2L) dat$Z[, 2L] <- 20 * dat$Z[, 2L]
      prepared <- pc_prepare_initial(dat$y, dat$X, dat$Z, dat$cluster_id, tau, diag(lambda, q),
                                      rep(0, 6L), fit_always_included = FALSE)
      conic <- pc_orthogonal_path(prepared, tau, 0L, gap_tolerance = 1e-9, solver = jp_joint_qp)
      reference <- pc_orthogonal_path(prepared, tau, 0L, gap_tolerance = 1e-9, max_iterations = 500000L)
      delta <- drop(dat$X %*% (conic$coefficients[, 1L] - reference$coefficients[, 1L]))
      distance <- sum(prepared$inverse_size * delta^2) / prepared$n
      bound <- 4 * prepared$n * (max(conic$fitted$duality_gap, 0) + max(reference$fitted$duality_gap, 0) + 1e-12)
      if (distance > bound || abs(conic$fitted$primal - reference$fitted$primal) > 2e-8) stop("The selected convex solvers disagree beyond their certificates.")
      records[[length(records) + 1L]] <- data.frame(q = q, tau = tau, lambda = lambda,
                                                   conic_gap = conic$fitted$duality_gap,
                                                   reference_gap = reference$fitted$duality_gap,
                                                   prediction_distance = distance, certificate_bound = bound)
    }
  }
}
write.csv(do.call(rbind, records), "results/scratch/joint_profile_manuscript_20261026/conic_checks.csv", row.names = FALSE)
cat("Passed", length(records), "selected convex solver comparisons\n")
