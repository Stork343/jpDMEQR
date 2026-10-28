library(data.table)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Provide a complete collected simulation file.")
source <- readRDS(args[1L])
metrics <- as.data.table(source$metrics)
coordinates <- as.data.table(source$coordinates)
keys <- c("experiment_id", "replicate", "method_id")
coordinates <- merge(coordinates, metrics[, c(keys, "optimizer_converged"), with = FALSE], by = keys, all.x = TRUE, sort = FALSE)
if (anyNA(coordinates$optimizer_converged)) stop("Coordinate records lack their fitting status.")
inference_methods <- unique(metrics$method_id[startsWith(metrics$method_id, "JOINT-PROFILE") |
  metrics$method_id %in% c("SQR-DEBIASED-IID", "BIAS-ADJ-LQMM")])
intervals <- coordinates[method_id %in% inference_methods]
intervals[, available := optimizer_converged & !is.na(feasible) & feasible &
  is.finite(beta_tilde) & is.finite(estimated_se) & estimated_se > 0 &
  is.finite(ci_lower) & is.finite(ci_upper) & is.finite(target_value)]
if (anyNA(intervals$available) || any(intervals$available &
  (is.na(intervals$covered) | is.na(intervals$covered_secondary)))) stop("An available interval lacks its recorded coverage.")
intervals[, covered_all := available & !is.na(covered) & covered]
intervals[, covered90_all := available & !is.na(covered_secondary) & covered_secondary]
intervals[, error := beta_tilde - target_value]
intervals[, rejected_zero := available & (ci_lower > 0 | ci_upper < 0)]
if (any(intervals$available & (intervals$covered != (intervals$ci_lower <= intervals$target_value & intervals$ci_upper >= intervals$target_value)), na.rm = TRUE))
  stop("Recorded interval coverage does not match the stored target and endpoints.")
mean_available <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
mcse_available <- function(x) if (sum(is.finite(x)) > 1L) sd(x[is.finite(x)]) / sqrt(sum(is.finite(x))) else NA_real_
binomial_bounds <- function(successes, attempted) {
  c(lower = if (successes == 0) 0 else qbeta(.025, successes, attempted - successes + 1),
    upper = if (successes == attempted) 1 else qbeta(.975, successes + 1, attempted - successes))
}
coordinate_summary <- intervals[, {
  valid <- available
  bounds <- binomial_bounds(sum(covered_all), .N)
  values <- error[valid]
  empirical_sd <- if (length(values) > 1L) sd(values) else NA_real_
  list(replications = .N, available_count = sum(valid), availability = mean(valid),
    coverage_all = mean(covered_all), coverage_mcse = sd(as.numeric(covered_all)) / sqrt(.N),
    coverage_lower = bounds[1L], coverage_upper = bounds[2L],
    coverage_available = if (any(valid)) mean(covered[valid]) else NA_real_,
    coverage90_all = mean(covered90_all), mean_error = mean_available(values),
    rejection_zero_all = mean(rejected_zero),
    mean_error_mcse = mcse_available(values), empirical_sd = empirical_sd,
    mean_se = mean_available(estimated_se[valid]),
    se_sd_ratio = mean_available(estimated_se[valid]) / empirical_sd,
    mean_length = mean_available(interval_length[valid]),
    length_mcse = mcse_available(interval_length[valid]),
    target_mcse = if (length(target_mc_se)) max(target_mc_se) else NA_real_)
}, by = .(experiment_id, method_id, coordinate, coordinate_type, target_type, target_value)]
replicate_groups <- intervals[, .(availability = mean(available), coverage = mean(covered_all),
  coverage90 = mean(covered90_all), length_total = mean(ifelse(available, interval_length, 0))),
  by = .(experiment_id, method_id, replicate, coordinate_type)]
group_summary <- replicate_groups[, {
  rate <- mean(availability)
  conditional <- if (is.finite(rate) && rate > 0) mean(coverage) / rate else NA_real_
  length_mean <- if (is.finite(rate) && rate > 0) mean(length_total) / rate else NA_real_
  list(replications = .N, availability = rate,
    coverage_all = mean(coverage), coverage_mcse = sd(coverage) / sqrt(.N),
    coverage_available = conditional,
    coverage_available_mcse = if (is.finite(rate) && rate > 0) sd(coverage - conditional * availability) / (sqrt(.N) * rate) else NA_real_,
    coverage90_all = mean(coverage90), coverage90_mcse = sd(coverage90) / sqrt(.N),
    mean_length = length_mean,
    length_mcse = if (is.finite(rate) && rate > 0) sd(length_total - length_mean * availability) / (sqrt(.N) * rate) else NA_real_)
},
  by = .(experiment_id, method_id, coordinate_type)]
predictive_summary <- metrics[, {
  valid <- optimizer_converged & is.finite(pinball_loss) & is.finite(quantile_calibration)
  losses <- pinball_loss[valid]
  calibration <- quantile_calibration[valid]
  list(replications = .N, converged_count = sum(optimizer_converged), predictive_count = sum(valid),
    availability = mean(valid), pinball_loss = mean_available(losses),
    pinball_loss_mcse = mcse_available(losses),
    calibration = mean_available(calibration), calibration_mcse = mcse_available(calibration),
    selected_size = mean_available(selected_size[valid]),
    runtime_median = median(runtime_sec), runtime_maximum = max(runtime_sec))
},
  by = .(experiment_id, method_id)]
estimation_summary <- metrics[, {
  valid <- optimizer_converged & is.finite(l2_error) & is.finite(l1_error)
  list(replications = .N, available_count = sum(valid),
    l2_error = mean_available(l2_error[valid]), l2_error_mcse = mcse_available(l2_error[valid]),
    l1_error = mean_available(l1_error[valid]), l1_error_mcse = mcse_available(l1_error[valid]),
    true_positive_rate = mean_available(tpr[valid]), false_discovery_proportion = mean_available(fdp[valid]),
    exact_support_probability = mean_available(as.numeric(exact_support[valid])))
}, by = .(experiment_id, method_id)]
output <- dirname(args[1L])
fwrite(coordinate_summary, file.path(output, "coordinate_summary.csv"))
fwrite(group_summary, file.path(output, "coordinate_group_summary.csv"))
fwrite(predictive_summary, file.path(output, "predictive_summary.csv"))
fwrite(estimation_summary, file.path(output, "estimation_summary.csv"))
saveRDS(replicate_groups, file.path(output, "replicate_groups.rds"))
cat("Summarised", nrow(coordinate_summary), "coordinate groups and", nrow(predictive_summary), "method settings\n")
