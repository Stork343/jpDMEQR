options(width = 180)
args <- commandArgs(trailingOnly = TRUE)
run_id <- if (length(args)) args[1L] else "paired_roots_independent"
root <- file.path("results/scratch/score_validity_20261026", run_id)
files <- sort(list.files(root, pattern = "^result[.]rds$", recursive = TRUE, full.names = TRUE))
if (!length(files)) stop("No completed evaluations.")
fits <- lapply(files, readRDS)
versions <- unique(vapply(fits, `[[`, character(1), "method_version"))
if (length(versions) != 1L || !versions %in% c("honest_paired_roots_20261026", "honest_profile_bootstrap_20261026", "honest_joint_profile_20261026", "joint_profile_design_20261026")) stop("Mixed or unknown method identifiers.")
table <- do.call(rbind, lapply(fits, `[[`, "table"))
table$type <- ifelse(table$truth == 0, "zero", "nonzero")
identity <- paste(table$case_id, table$replicate, table$coordinate)
if (anyDuplicated(identity)) stop("Duplicate evaluation coordinates.")
coordinate <- lapply(split(table, paste(table$case_id, table$coordinate)), function(data) {
  coverage <- mean(data$covered)
  ci <- binom.test(sum(data$covered), nrow(data))$conf.int
  data.frame(case_id = data$case_id[1L], coordinate = data$coordinate[1L], type = data$type[1L],
             replications = nrow(data), coverage = coverage, mcse = sqrt(coverage * (1 - coverage) / nrow(data)),
             mc_lower = ci[1L], mc_upper = ci[2L],
             mean_error = mean(data$beta_tilde - data$truth), empirical_sd = sd(data$beta_tilde),
             mean_se = mean(data$se), se_sd_ratio = mean(data$se) / sd(data$beta_tilde),
             interval_length = mean(data$upper - data$lower), rejection = mean(data$reject_zero),
             z_mean = mean(data$z_truth), z_sd = sd(data$z_truth))
})
coordinate <- do.call(rbind, coordinate)
aggregate_coverage <- function(data) {
  grouped <- aggregate(covered ~ case_id + replicate, data, mean)
  strata <- split(grouped, grouped$case_id)
  total <- nrow(grouped)
  variance <- sum(vapply(strata, function(group) {
    if (nrow(group) < 2L) return(NA_real_)
    (nrow(group) / total)^2 * var(group$covered) / nrow(group)
  }, numeric(1)))
  coverage <- mean(grouped$covered)
  mcse <- sqrt(variance)
  data.frame(datasets = total, intervals = nrow(data), coverage = coverage, mcse = mcse,
             mc_lower = max(0, coverage - qnorm(0.975) * mcse),
             mc_upper = min(1, coverage + qnorm(0.975) * mcse),
             z_mean = mean(data$z_truth), z_sd = sd(data$z_truth),
             rejection = mean(data$reject_zero), interval_length = mean(data$upper - data$lower))
}
overall <- do.call(rbind, lapply(split(table, table$type), aggregate_coverage))
overall$type <- rownames(overall)
scenario <- do.call(rbind, lapply(split(table, paste(table$case_id, table$type)), function(data) {
  cbind(data.frame(case_id = data$case_id[1L], type = data$type[1L]), aggregate_coverage(data))
}))
timings <- do.call(rbind, lapply(fits, function(fit) {
  data.frame(case_id = fit$request$config$case_id, replicate = fit$request$replicate,
             training = unname(fit$timing["training"]), inference = unname(fit$timing["inference"]),
             total = unname(fit$timing["total"]),
             max_gap = if (fit$method_version == "honest_paired_roots_20261026") {
               max(vapply(fit$training, function(fold) max(vapply(fold$refits, `[[`, numeric(1), "gap")), numeric(1)))
             } else NA_real_,
             max_kkt = if (fit$method_version == "honest_profile_bootstrap_20261026") {
               max(vapply(fit$details, function(detail) max(detail$bootstrap[, "max_kkt"]), numeric(1)))
             } else if (fit$method_version %in% c("honest_joint_profile_20261026", "joint_profile_design_20261026")) max(fit$table$max_kkt) else NA_real_,
             largest_support = max(vapply(fit$training, function(fold) max(lengths(lapply(fold$graphs, `[[`, "support"))), integer(1))))
}))
failures <- list.files(root, pattern = "^failure[.]rds$", recursive = TRUE, full.names = TRUE)
failure_table <- data.frame()
if (length(failures)) {
  failure_table <- do.call(rbind, lapply(failures, function(path) {
    failure <- readRDS(path)
    request <- failure$request
    resolved <- any(table$case_id == request$config$case_id & table$replicate == request$replicate)
    data.frame(case_id = request$config$case_id, replicate = request$replicate,
               error = conditionMessage(failure$condition), resolved = resolved, file = path)
  }))
}
directory <- file.path("results/tables/score_validity_20261026", run_id)
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
write.csv(table, file.path(directory, "coordinates.csv"), row.names = FALSE)
write.csv(coordinate, file.path(directory, "coordinate_summary.csv"), row.names = FALSE)
write.csv(overall, file.path(directory, "overall.csv"), row.names = FALSE)
write.csv(scenario, file.path(directory, "scenario_summary.csv"), row.names = FALSE)
write.csv(timings, file.path(directory, "timings.csv"), row.names = FALSE)
write.csv(failure_table, file.path(directory, "failures.csv"), row.names = FALSE)
cat("Completed datasets", length(fits), "intervals", nrow(table), "failure records", length(failures), "\n")
print(overall, row.names = FALSE)
print(scenario[, c("case_id", "type", "datasets", "coverage", "z_mean", "z_sd")], row.names = FALSE)
