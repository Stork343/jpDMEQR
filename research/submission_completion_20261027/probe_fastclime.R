source("research/submission_completion_20261027/load.R")
invisible(loadNamespace("fastclime"))
output <- "results/scratch/submission_completion_20261027/fastclime_probe"
if (dir.exists(output)) stop("The fastclime numerical probe already exists.")
dir.create(output, recursive = TRUE)
reference <- readRDS("results/scratch/submission_comparator_probes_20261027/SQR-DEBIASED-IID_A01_matrix_free_fista/result.rds")$fit
H <- reference$fit_object$hessian
p <- nrow(H)
records <- list()
for (k in 1:5) {
  canonical <- numeric(p)
  canonical[k] <- 1
  steps <- 4L * p
  started <- proc.time()[[3L]]
  path <- .C("dantzig", as.double(H), as.double(canonical), beta = numeric(p * steps),
    as.integer(p), as.double(.05), as.integer(steps), gamma = numeric(steps), PACKAGE = "fastclime")
  elapsed <- proc.time()[[3L]] - started
  valid <- which(path$gamma > 0)
  gamma_path <- path$gamma[valid]
  coefficients <- matrix(path$beta, p, steps)[, valid, drop = FALSE]
  saveRDS(list(gamma = gamma_path, beta = coefficients, elapsed = elapsed),
    file.path(output, paste0("coordinate_", k, ".rds")))
  for (gamma in c(.15, .3, .6)) {
    if (min(gamma_path) > gamma) stop("The fastclime path did not reach the requested penalty.")
    upper <- max(which(gamma_path >= gamma))
    lower <- upper + 1L
    weight <- (gamma - gamma_path[lower]) / (gamma_path[upper] - gamma_path[lower])
    omega <- weight * coefficients[, upper] + (1 - weight) * coefficients[, lower]
    residual <- max(abs(H %*% omega - canonical))
    records[[length(records) + 1L]] <- data.frame(coordinate = k, gamma = gamma, elapsed = elapsed,
      path_length = length(valid), endpoint = min(gamma_path), residual = residual, l1_norm = sum(abs(omega)))
  }
  write.csv(do.call(rbind, records), file.path(output, "checks.csv"), row.names = FALSE)
}
print(do.call(rbind, records))
