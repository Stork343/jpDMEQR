source("research/score_validity_20261026/load_joint.R")
source("research/joint_profile_manuscript_20261026/fixed_inference.R")
records <- list()
for (q in 1:2) {
  for (tau in c(0.25, 0.5, 0.75)) {
    dat <- generate_profile_qr_data_v2(n_clusters = 80L, p = 8L, s = 4L, tau = tau,
                                       q = q, seed = 202610260L + q * 10L + as.integer(4 * tau))
    Lambda <- if (q == 1L) diag(1L) else matrix(c(1, 0.2, 0.2, 0.8), 2L)
    fit <- jp_fit_fixed(dat$y, dat$X, dat$Z, dat$cluster_id, tau, c(1L, 2L, 6L), Lambda)
    for (j in 1:3) {
      k <- c(1L, 2L, 6L)[j]
      order <- c(k, setdiff(1:8, k))
      data <- sb_data(dat$y, dat$X[, order, drop = FALSE], dat$Z, dat$cluster_id)
      reference <- ds_pair_small(data, numeric(8L), tau, 80^(-0.3), Lambda)
      se <- sqrt(mean((reference$influence - mean(reference$influence))^2) / 80)
      estimate_error <- abs(reference$estimate - fit$table$beta_hat[j])
      se_error <- abs(se - fit$table$se[j])
      if (max(estimate_error, se_error) > 2e-6) stop("The fixed design interface differs from its coordinate reference.")
      records[[length(records) + 1L]] <- data.frame(q = q, tau = tau, coordinate = k,
                                                   estimate_error = estimate_error, se_error = se_error)
    }
    if (min(eigen(fit$covariance, symmetric = TRUE, only.values = TRUE)$values) < -1e-12) stop("Invalid joint covariance.")
  }
}
dir.create("results/scratch/joint_profile_manuscript_20261026", recursive = TRUE, showWarnings = FALSE)
write.csv(do.call(rbind, records), "results/scratch/joint_profile_manuscript_20261026/fixed_checks.csv", row.names = FALSE)
cat("Passed", length(records), "fixed design coordinate comparisons\n")
