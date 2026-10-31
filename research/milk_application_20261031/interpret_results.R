library(data.table)
library(yaml)
library(quantreg)
config <- read_yaml("config/application_milk_20261031.yml")
root <- file.path(config$results_directory, config$execution_version)
inference <- fread("results/preflight/milk_scs_20261031/inference.csv")
reference <- inference[setting == "reference"][order(tau)]
if (!identical(reference$tau, c(0.25, 0.5, 0.75))) stop("The registered reference quantiles are incomplete.")
influence <- matrix(NA_real_, reference$n[1L], nrow(reference))
subjects <- NULL
fit_commits <- character(nrow(reference))
for (j in seq_len(nrow(reference))) {
  request <- readRDS(file.path(root, reference$case_id[j], "request.rds"))
  fit_commits[j] <- request$source_commit
  fit <- readRDS(file.path(root, reference$case_id[j], "fit.rds"))
  current_subjects <- unlist(lapply(fit$prepared$folds, `[[`, "validation_subjects"), use.names = FALSE)
  if (anyDuplicated(current_subjects) || length(current_subjects) != reference$n[j])
    stop("The fitted cow identifiers do not define one contribution per cow.")
  if (is.null(subjects)) subjects <- current_subjects
  index <- match(subjects, current_subjects)
  if (anyNA(index)) stop("The reference estimates use different cows.")
  influence[, j] <- fit$details[[1L]]$influence[index]
}
n <- length(subjects)
covariance <- crossprod(influence) / n^2
if (max(abs(sqrt(diag(covariance)) - reference$se)) > 1e-12)
  stop("The joint covariance does not reproduce the saved standard errors.")
contrast_matrix <- rbind(c(-1, 1, 0), c(-1, 0, 1), c(0, -1, 1))
rownames(contrast_matrix) <- c("q50_minus_q25", "q75_minus_q25", "q75_minus_q50")
estimate <- drop(contrast_matrix %*% reference$beta_tilde)
contrast_covariance <- contrast_matrix %*% covariance %*% t(contrast_matrix)
se <- sqrt(diag(contrast_covariance))
direct_variance <- colSums((influence %*% t(contrast_matrix))^2) / n^2
if (max(abs(se^2 - direct_variance)) > 1e-12 || any(!is.finite(se)) || any(se <= 0))
  stop("The contrast covariance identity failed.")
critical <- qnorm(1 - 0.05 / (2 * nrow(contrast_matrix)))
contrasts <- data.table(contrast = rownames(contrast_matrix), estimate = estimate, se = se,
  simultaneous_lower = estimate - critical * se, simultaneous_upper = estimate + critical * se,
  p_value = 2 * pnorm(-abs(estimate / se)))
contrasts[, bonferroni_p := p.adjust(p_value, method = "bonferroni")]
global_covariance <- contrast_covariance[1:2, 1:2]
if (min(eigen(global_covariance, symmetric = TRUE, only.values = TRUE)$values) <= 0)
  stop("The global contrast covariance is not identified.")
statistic <- sum(estimate[1:2] * solve(global_covariance, estimate[1:2]))
global_test <- data.table(hypothesis = "theta_q25_equals_theta_q50_equals_theta_q75",
  statistic = statistic, degrees_freedom = 2L, p_value = pchisq(statistic, 2L, lower.tail = FALSE))
cohort <- readRDS(file.path(config$results_directory, "cohorts.rds"))$eligible
observations <- data.table(cluster_id = cohort$cluster_id, scs = cohort$y,
  A2_milk = cohort$X[, "A2_milk"])
observations[, visits := .N, by = cluster_id]
observations[, weight := 1 / visits]
groups <- list(all = seq_len(nrow(observations)), A1 = which(observations$A2_milk == 0),
  A2 = which(observations$A2_milk == 1), single_control = which(observations$visits == 1L),
  two_controls = which(observations$visits == 2L))
descriptive <- rbindlist(lapply(names(groups), function(group) {
  current <- observations[groups[[group]]]
  quantiles <- as.numeric(coef(quantreg::rq(scs ~ 1, weights = weight, data = current,
    tau = c(0.25, 0.5, 0.75), method = "fn")))
  if (length(quantiles) != 3L || any(!is.finite(quantiles)) || any(diff(quantiles) < 0))
    stop("The observed weighted quantiles are invalid.")
  cows <- unique(current[, .(cluster_id, A2_milk)])
  data.table(group = group, cows = nrow(cows), samples = nrow(current),
    A2_fraction = mean(cows$A2_milk), cow_weighted_mean = weighted.mean(current$scs, current$weight),
    cow_weighted_q25 = quantiles[1L], cow_weighted_q50 = quantiles[2L], cow_weighted_q75 = quantiles[3L])
}))
sensitivity <- inference[, .(estimate_min = min(beta_tilde), estimate_max = max(beta_tilde),
  all_simultaneous_intervals_include_zero = all(simultaneous_lower <= 0 & simultaneous_upper >= 0)), by = setting]
bandwidth <- inference[setting %in% c("bandwidth_075", "bandwidth_125")]
bandwidth[, reference_estimate := reference$beta_tilde[match(tau, reference$tau)]]
paired_median <- inference[setting == "paired" & tau == 0.5]
interpretation <- list(maximum_bandwidth_change = max(abs(bandwidth$beta_tilde - bandwidth$reference_estimate)),
  paired_median_estimate = paired_median$beta_tilde,
  paired_median_lower = paired_median$simultaneous_lower,
  paired_median_upper = paired_median$simultaneous_upper,
  all_registered_intervals_include_zero = all(inference$simultaneous_lower <= 0 & inference$simultaneous_upper >= 0),
  single_control_mean = descriptive[group == "single_control", cow_weighted_mean],
  two_control_mean = descriptive[group == "two_controls", cow_weighted_mean],
  marginal_A1_median = descriptive[group == "A1", cow_weighted_q50],
  marginal_A2_median = descriptive[group == "A2", cow_weighted_q50])
output <- "results/preflight/milk_scs_20261031/interpretation"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
fwrite(contrasts, file.path(output, "quantile_contrasts.csv"))
fwrite(global_test, file.path(output, "quantile_equality_test.csv"))
fwrite(descriptive, file.path(output, "observed_distribution.csv"))
fwrite(sensitivity, file.path(output, "sensitivity_interpretation.csv"))
fwrite(as.data.table(covariance), file.path(output, "reference_covariance.csv"))
jsonlite::write_json(interpretation, file.path(output, "interpretation_values.json"),
  auto_unbox = TRUE, pretty = TRUE, digits = NA)
claims <- data.table(question = c("adjusted_milk_type_association", "variation_across_quantiles",
  "bandwidth_sensitivity", "two_control_population", "computational_feasibility"),
  source = c("../inference.csv", "quantile_contrasts.csv;quantile_equality_test.csv",
    "interpretation_values.json", "observed_distribution.csv;../inference.csv", "../completion_manifest.json"),
  conclusion = c("三个同时区间均包括零，关联方向尚未确定", "整体检验与全部成对区间未检出分位数差异",
    "报告带宽变化引起的最大系数变化", "同时报告估计中心变化和奶牛样本的结局分布",
    "实际数据上完成完整高维坐标推断并报告运行时间"))
fwrite(claims, file.path(output, "claim_mapping.csv"))
jsonlite::write_json(list(analysis_role = "supplementary_interpretation", new_high_dimensional_fits = 0L,
  source_cases = reference$case_id, coordinate_interval_family = 3L, contrast_interval_family = 3L,
  covariance_standard_error_error = max(abs(sqrt(diag(covariance)) - reference$se)),
  contrast_variance_identity_error = max(abs(se^2 - direct_variance)),
  source_fit_commits = unique(fit_commits),
  source_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  computed = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE)),
  file.path(output, "manifest.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
print(contrasts)
print(global_test)
print(descriptive)
