source("research/score_validity_20261026/load_joint.R")
records <- list()
record <- function(label, error, tolerance) {
  records[[length(records) + 1L]] <<- data.frame(check = label, error = error, tolerance = tolerance)
  if (!is.finite(error) || error > tolerance) stop(paste("Failed joint inference check", label, error))
}
for (q in 1:2) {
  for (tau in c(0.25, 0.5, 0.75)) {
    dat <- generate_profile_qr_data_v2(n_clusters = 120L, p = 100L, s = 5L, tau = tau, q = q,
                                       seed = 1200600000L + 10L * q + as.integer(4 * tau))
    fit <- ds_fit_profile(dat$y, dat$X, dat$Z, dat$cluster_id, tau, c(1L, 6L), keep_prepared = TRUE)
    for (f in 1:2) {
      training <- fit$training[[f]]
      record(paste0("subjects_", q, "_", tau, "_", f),
              length(intersect(training$training_subjects, training$validation_subjects)), 0)
    }
    for (j in 1:2) {
      moment <- 0
      for (f in 1:2) {
        small <- fit$details[[j]]$fits[[f]]
        full <- small$full$projected_scores - mean(small$full$projected_scores)
        half <- small$half$projected_scores - mean(small$half$projected_scores)
        moment <- moment + 4 * sum(half^2) + sum(full^2) - 4 * sum(full * half)
      }
      record(paste0("paired_covariance_", q, "_", tau, "_", j),
              abs(fit$table$se[j]^2 - moment / fit$prepared$n^2), 1e-12)
      record(paste0("fitted_kkt_", q, "_", tau, "_", j), fit$table$max_kkt[j], 1e-8)
    }
    scaled <- dat$X
    scaled[, 1L] <- 3 * scaled[, 1L]
    scaled[, 4L] <- 0.5 * scaled[, 4L]
    other <- ds_fit_profile(dat$y, scaled, dat$Z, dat$cluster_id, tau, c(1L, 6L))
    target_scale <- c(3, 1)
    for (column in c("beta_tilde", "se", "lower", "upper")) {
      record(paste0("units_", q, "_", tau, "_", column),
              max(abs(other$table[[column]] * target_scale - fit$table[[column]])), 2e-5)
    }
    reordered <- rev(seq_along(dat$y))
    ordered <- ds_fit_profile(dat$y[reordered], dat$X[reordered, , drop = FALSE],
                               dat$Z[reordered, , drop = FALSE], dat$cluster_id[reordered], tau, c(1L, 6L))
    record(paste0("row_order_", q, "_", tau),
            max(abs(as.matrix(ordered$table[, 2:5]) - as.matrix(fit$table[, 2:5]))), 2e-5)
  }
}
result <- do.call(rbind, records)
write.csv(result, "results/scratch/score_validity_20261026/design_joint_checks.csv", row.names = FALSE)
cat("Passed", nrow(result), "joint inference checks\n")
