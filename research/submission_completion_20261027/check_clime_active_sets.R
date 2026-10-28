source("research/submission_completion_20261027/load.R")
source("research/submission_completion_20261027/clime_flare.R")
source("research/submission_completion_20261027/clime_certified.R")
H <- readRDS("results/scratch/submission_completion_20261027/precision_row_checks/H.rds")
gamma <- .152485166387
started <- proc.time()[[3L]]
warm <- submission_flare_matrix(H, gamma, max_iterations = 1000L)
elapsed <- proc.time()[[3L]] - started
checks <- lapply(seq_len(nrow(H)), function(k) {
  candidate <- submission_clime_certificate(H, k, gamma, warm$matrices[[1L]][, k])
  data.frame(coordinate = k, verified = candidate$verified,
    support = sum(abs(warm$matrices[[1L]][, k]) > 1e-6),
    residual = candidate$residual %||% NA_real_, dual_residual = candidate$dual_residual %||% NA_real_,
    dual_gap = candidate$dual_gap %||% NA_real_)
})
records <- do.call(rbind, checks)
saveRDS(list(H = H, gamma = gamma, warm = warm, elapsed = elapsed, checks = records),
  "results/scratch/submission_completion_20261027/clime_active_set_check.rds")
write.csv(records, "results/scratch/submission_completion_20261027/clime_active_set_checks.csv", row.names = FALSE)
cat("Warm direction seconds", elapsed, "certified", sum(records$verified), "of", nrow(records), "\n")
print(summary(records))
