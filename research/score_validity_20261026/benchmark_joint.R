source("research/score_validity_20261026/load_joint.R")
root <- "results/scratch/score_validity_20261026/design_scaling"
if (dir.exists(root)) stop("The scaling run already exists.")
dir.create(root, recursive = TRUE)
saveRDS(list(started = Sys.time(), session = sessionInfo(), system = Sys.info(),
             workers = 1L, replications = 20L, coordinates = c(1L, 2L, 6L, 7L),
             method_version = "joint_profile_design_20261026"), file.path(root, "request.rds"))
records <- list()
for (p in c(500L, 1000L, 2000L, 10000L)) {
  for (replicate in seq_len(20L)) {
    dat <- generate_profile_qr_data_v2(n_clusters = 200L, p = p, s = 5L, tau = 0.5, q = 1L,
                                       seed = 1900000000L + 1000L * p + replicate)
    started <- proc.time()
    fit <- ds_fit_profile(dat$y, dat$X, dat$Z, dat$cluster_id, 0.5, c(1L, 2L, 6L, 7L))
    used <- proc.time() - started
    records[[length(records) + 1L]] <- data.frame(p = p, replicate = replicate,
                                                 elapsed = used["elapsed"], user = used["user.self"],
                                                 system = used["sys.self"], training = fit$timing["training"],
                                                 inference = fit$timing["inference"], max_kkt = max(fit$table$max_kkt),
                                                 largest_support = max(vapply(fit$training, function(fold) max(lengths(lapply(fold$graphs, `[[`, "support"))), integer(1))))
    write.csv(do.call(rbind, records), file.path(root, "timings.csv"), row.names = FALSE)
  }
  cat("Completed dimension", p, "\n")
}
timings <- do.call(rbind, records)
summary <- do.call(rbind, lapply(split(timings, timings$p), function(data) {
  data.frame(p = data$p[1L], replications = nrow(data), median = median(data$elapsed),
             p90 = quantile(data$elapsed, .9, names = FALSE), maximum = max(data$elapsed),
             median_training = median(data$training), median_inference = median(data$inference),
             max_kkt = max(data$max_kkt), maximum_selected_size = max(data$largest_support))
}))
write.csv(summary, file.path(root, "summary.csv"), row.names = FALSE)
saveRDS(list(finished = Sys.time(), completed = nrow(timings)), file.path(root, "complete.rds"))
print(summary)
