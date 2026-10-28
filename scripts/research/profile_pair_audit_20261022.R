args <- base::commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1]] else ".", mustWork = TRUE)
out <- file.path(root, "results", "scratch", "profile_score_audit_20261022")
if (identical(Sys.getenv("JPDMEQR_REUSE_RESEARCH_RESULTS"), "TRUE")) {
  records <- utils::read.csv(file.path(out, "paired_components.csv"))
  stopifnot(nrow(records) == 18L,
            max(records$score_fd_error) < 2e-4,
            max(records$hessian_fd_error) < 3e-3,
            max(records$covariance_identity_error) < 1e-12)
  print(vapply(records[c("score_fd_error", "hessian_fd_error",
                         "covariance_identity_error", "inverse_defect")], max, numeric(1)))
  quit(save = "no", status = 0L)
}
stage <- "initialisation"
solver <- Sys.getenv("JPDMEQR_CLUSTER_SOLVER", unset = "bfgs")
started <- as.character(Sys.time())
attempt_path <- file.path(out, paste0("pair_attempt_", solver, ".csv"))
record_attempt <- function(status) {
  utils::write.csv(data.frame(solver = solver, started = started,
                             finished = as.character(Sys.time()),
                             stage = stage, status = status),
                   attempt_path, row.names = FALSE)
}
options(error = function() {
  record_attempt("failed")
  writeLines(c(paste("stage:", stage), geterrmessage()),
             file.path(out, paste0("PAIR_FAILED_", solver, ".txt")))
  quit(save = "no", status = 1L)
})
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2")
records <- list()
index <- 0L
for (tau in c(0.25, 0.5, 0.75)) {
  dat <- readRDS(file.path(out, paste0("cluster_reference_tau", tau, ".rds")))$data
  for (h in c(0.2, 0.1, 0.05)) {
    stage <- paste("tau", tau, "h", h)
    evaluate <- function(beta) profile_bandwidth_pair_v2(
      dat$y, dat$X, dat$Z, dat$cluster_id, beta, tau, h,
      identity_check = TRUE
    )
    beta <- as.numeric(dat$beta0)
    pair <- evaluate(beta)
    omega <- solve(pair$hessian)
    interval <- profile_bandwidth_pair_wald_v2(pair, beta, omega, 1:2)
    delta <- 2e-4
    for (coordinate in 1:2) {
      direction <- replace(numeric(2), coordinate, delta)
      plus <- evaluate(beta + direction)
      minus <- evaluate(beta - direction)
      Hfd <- (plus$score - minus$score) / (2 * delta)
      gfd <- (plus$objective - minus$objective) / (2 * delta)
      half <- as.numeric(pair$components_half_h$cluster_scores %*% omega[, coordinate])
      full <- as.numeric(pair$components_h$cluster_scores %*% omega[, coordinate])
      half <- half - mean(half)
      full <- full - mean(full)
      expanded <- 4 * mean(half^2) + mean(full^2) - 4 * mean(half * full)
      index <- index + 1L
      records[[index]] <- data.frame(
        tau = tau, h = h, coordinate = coordinate,
        score_fd_error = abs(gfd - pair$score[coordinate]),
        hessian_fd_error = max(abs(Hfd - pair$hessian[, coordinate])),
        variance = interval$se[coordinate]^2 * pair$n_clusters,
        expanded_variance = expanded,
        cross_bandwidth_covariance = mean(half * full),
        covariance_identity_error = abs(interval$se[coordinate]^2 * pair$n_clusters - expanded),
        inverse_defect = interval$inverse_defect[coordinate],
        hessian_min_eigenvalue = pair$hessian_min_eigenvalue,
        nuisance_gradient_max = pair$max_nuisance_gradient
      )
    }
    utils::write.csv(do.call(rbind, records),
                     file.path(out, paste0("paired_components_", solver, ".csv")),
                     row.names = FALSE)
  }
}
records <- do.call(rbind, records)
stopifnot(max(records$score_fd_error) < 2e-4,
          max(records$hessian_fd_error) < 3e-3,
          max(records$covariance_identity_error) < 1e-12)
utils::write.csv(records, file.path(out, "paired_components.csv"), row.names = FALSE)
record_attempt("complete")
print(vapply(records[c("score_fd_error", "hessian_fd_error",
                       "covariance_identity_error", "inverse_defect")], max, numeric(1)))
