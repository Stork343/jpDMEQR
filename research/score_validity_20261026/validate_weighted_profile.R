source("research/score_validity_20261026/load_joint.R")
checks <- list()
record <- function(label, error, tolerance) {
  checks[[length(checks) + 1L]] <<- data.frame(check = label, error = error, tolerance = tolerance)
  if (!is.finite(error) || error > tolerance) stop(paste("Failed", label, error))
}
for (seed in 1:20) {
  for (q in 1:2) {
    for (tau in c(0.25, 0.5, 0.75)) {
      dat <- generate_profile_qr_data_v2(n_clusters = 40L, p = 5L, s = 3L, tau = tau, q = q,
                                         seed = 1200400000L + seed * 100L + q * 10L + as.integer(4 * tau))
      data <- sb_data(dat$y, dat$X, dat$Z, dat$cluster_id)
      Lambda <- if (q == 1L) diag(1) else matrix(c(1, 0.2, 0.2, 0.8), 2L)
      weights <- as.numeric(rmultinom(1L, size = data$n, prob = rep(1 / data$n, data$n)))
      fit <- ds_fit_small(data, dat$beta0, tau, 0.4, Lambda, weights)
      gram <- crossprod(data$X, data$X * data$weight) / data$n
      indices <- unlist(lapply(seq_len(data$n), function(i) rep(data$rows[[i]], weights[i])), use.names = FALSE)
      ids <- unlist(lapply(seq_len(data$n), function(i) {
        unlist(lapply(seq_len(weights[i]), function(j) rep(paste0(i, "_", j), length(data$rows[[i]]))), use.names = FALSE)
      }), use.names = FALSE)
      expanded <- sb_data(data$y[indices], data$X[indices, , drop = FALSE],
                          data$Z[indices, , drop = FALSE], ids)
      duplicate <- ds_fit_small(expanded, dat$beta0, tau, 0.4, Lambda, gram = gram)
      label <- paste(seed, q, tau, sep = "_")
      record(paste0("weighted_multiplicity_", label), max(abs(fit$beta - duplicate$beta)), 2e-6)
      profile <- ds_profile(data, fit$beta, tau, 0.4, Lambda)
      F <- sb_geometry_cpp(data$X, data$Z, profile$phi, data$order, data$pointers, Lambda, TRUE, weights)
      H <- sb_geometry_cpp(data$X, data$Z, profile$phi, data$order, data$pointers, Lambda, FALSE, weights)
      record(paste0("weighted_factor_", label), max(abs(crossprod(F) / data$n - H)), 1e-10)
      gradient <- function(beta) {
        prof <- ds_profile(data, beta, tau, 0.4, Lambda)
        -drop(crossprod(data$X, prof$psi * data$weight * weights[data$group])) / data$n
      }
      direction <- c(1, -0.5, 0.2, 0.4, -0.1)
      derivative <- (gradient(fit$beta + 1e-5 * direction) - gradient(fit$beta - 1e-5 * direction)) / 2e-5
      record(paste0("weighted_derivative_", label), max(abs(derivative - H %*% direction)), 2e-5)
      record(paste0("newton_kkt_", label), fit$kkt, 1e-8)
    }
  }
}
write.csv(do.call(rbind, checks), "results/scratch/score_validity_20261026/weighted_profile_checks.csv", row.names = FALSE)
cat("Passed", length(checks), "weighted profile checks\n")
