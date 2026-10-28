source("research/score_validity_20261026/load_joint.R")
paths <- c(
  "results/scratch/score_validity_20261026/design_independent/N200_Q1_T25/replicate_001/result.rds",
  "results/scratch/score_validity_20261026/design_independent/N400_P2000/replicate_060/result.rds")
records <- list()
for (path in paths) {
  stored <- readRDS(path)
  cfg <- stored$request$config
  dat <- generate_profile_qr_data_v2(n_clusters = cfg$n_clusters, p = cfg$p, s = cfg$s,
                                     tau = cfg$tau, q = cfg$q, error_dist = cfg$error_dist,
                                     signal = cfg$signal, rho_x = cfg$rho_x, seed = stored$request$seed)
  repeated <- ds_fit_profile(dat$y, dat$X, dat$Z, dat$cluster_id, cfg$tau, c(1L, 2L, 6L, 7L))
  columns <- c("beta_tilde", "se", "lower", "upper")
  error <- max(abs(as.matrix(stored$table[, columns]) - as.matrix(repeated$table[, columns])))
  if (!is.finite(error) || error > 1e-8) stop("Reinstalled package does not reproduce the stored estimator.")
  records[[length(records) + 1L]] <- data.frame(case_id = cfg$case_id, seed = stored$request$seed, maximum_difference = error)
}
result <- do.call(rbind, records)
write.csv(result, "results/scratch/score_validity_20261026/reproduction_checks.csv", row.names = FALSE)
print(result)
