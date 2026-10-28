library(data.table)
output <- "results/tables/submission_20261028/collections/statistical_study"
source <- readRDS(file.path(output, "records.rds"))
metrics <- as.data.table(source$metrics)
groups <- readRDS(file.path(output, "replicate_groups.rds"))
core <- groups[method_id %in% c("JOINT-PROFILE", "JOINT-PROFILE-SPLIT")]
others <- groups[!method_id %in% c("JOINT-PROFILE", "JOINT-PROFILE-SPLIT")]
join <- c("experiment_id", "replicate", "coordinate_type")
paired <- merge(others, core, by = join, all.x = TRUE, suffixes = c("_comparison", "_joint"))
if (nrow(paired) != nrow(others) || anyNA(paired$method_id_joint)) stop("An inferential comparator lacks its paired joint estimate.")
coverage <- paired[, {
  difference <- coverage_joint - coverage_comparison
  list(replications = .N, difference = mean(difference),
    difference_mcse = sd(difference) / sqrt(.N),
    joint_availability = mean(availability_joint), comparison_availability = mean(availability_comparison))
}, by = .(experiment_id, method_id_comparison, coordinate_type)]
core_metrics <- metrics[method_id %in% c("JOINT-PROFILE", "JOINT-PROFILE-SPLIT"),
  .(experiment_id, replicate, joint_loss = pinball_loss, joint_calibration = quantile_calibration,
    joint_converged = optimizer_converged)]
other_metrics <- metrics[!method_id %in% c("JOINT-PROFILE", "JOINT-PROFILE-SPLIT")]
prediction <- merge(other_metrics, core_metrics, by = c("experiment_id", "replicate"), all.x = TRUE)
if (nrow(prediction) != nrow(other_metrics) || anyNA(prediction$joint_converged)) stop("An estimator lacks its paired joint record.")
prediction <- prediction[, {
  valid <- optimizer_converged & joint_converged & is.finite(pinball_loss) & is.finite(joint_loss)
  difference <- joint_loss[valid] - pinball_loss[valid]
  list(replications = .N, paired_available = sum(valid),
    difference = if (length(difference)) mean(difference) else NA_real_,
    difference_mcse = if (length(difference) > 1L) sd(difference) / sqrt(length(difference)) else NA_real_)
}, by = .(experiment_id, method_id)]
fwrite(coverage, file.path(output, "paired_coverage_differences.csv"))
fwrite(prediction, file.path(output, "paired_prediction_differences.csv"))
cat("Calculated", nrow(coverage), "paired coverage groups and", nrow(prediction), "prediction comparisons\n")
