library(data.table)
library(ggplot2)
root <- "results/raw/simulation/submission_20261028_gaussian_diagnostics"
request <- readRDS("results/raw/simulation/submission_20261027_main_joint/request.rds")
registry <- as.data.table(request$registry)[module == "A_theory_scaling"]
jobs <- as.data.table(request$jobs)[experiment_id %in% registry$experiment_id]
coordinates <- folds <- accounting <- vector("list", nrow(jobs))
for (i in seq_len(nrow(jobs))) {
  job <- jobs[i]
  path <- file.path(root, job$experiment_id, sprintf("replicate_%05d.rds", job$replicate))
  if (!file.exists(path)) stop("The Gaussian diagnostic study is incomplete.")
  result <- readRDS(path)
  if (result$experiment_id != job$experiment_id || result$replicate != job$replicate) stop("A diagnostic identity differs.")
  accounting[[i]] <- data.table(experiment_id = job$experiment_id, replicate = job$replicate, available = result$available)
  if (!result$available) next
  if (!identical(result$curvature$definition, "Gaussian_Stein_expected_innovation_times_profile_dual") ||
      nrow(result$coordinates) != 4L || nrow(result$folds) != 16L) stop("A Gaussian diagnostic definition differs.")
  coordinates[[i]] <- cbind(job[, .(experiment_id, replicate)], as.data.table(result$coordinates))
  folds[[i]] <- cbind(job[, .(experiment_id, replicate)], as.data.table(result$folds))
}
accounting <- rbindlist(accounting)
coordinates <- rbindlist(coordinates)
folds <- rbindlist(folds)
if (anyDuplicated(coordinates[, .(experiment_id, replicate, coordinate)]) ||
    anyDuplicated(folds[, .(experiment_id, replicate, coordinate, fold, bandwidth)])) stop("Duplicated theory diagnostics.")
coordinate_replications <- coordinates[, .(squared_remainder = mean(scaled_remainder^2),
  squared_relative_remainder = mean(remainder_over_se^2), maximum_target_shift = max(maximum_target_shift)),
  by = .(experiment_id, replicate, coordinate_type)]
coordinate_summary <- coordinate_replications[, {
  rms <- sqrt(mean(squared_remainder))
  relative <- sqrt(mean(squared_relative_remainder))
  list(available_replications = .N, remainder_rms = rms,
    remainder_rms_mcse = sd(squared_remainder) / (2 * sqrt(.N) * rms),
    relative_remainder_rms = relative, relative_remainder_mcse = sd(squared_relative_remainder) / (2 * sqrt(.N) * relative),
    maximum_target_shift = max(maximum_target_shift))
}, by = .(experiment_id, coordinate_type)]
fold_replications <- folds[, .(curvature_max_error = mean(curvature_max_error),
  curvature_operator_error = mean(curvature_operator_error),
  inverse_row_l2_error = mean(inverse_row_l2_error), coefficient_l2_error = mean(coefficient_l2_error),
  selected_dimension = mean(selected_dimension), maximum_kkt = max(kkt)), by = .(experiment_id, replicate)]
fold_summary <- fold_replications[, lapply(.SD, mean), by = experiment_id, .SDcols =
  c("curvature_max_error", "curvature_operator_error", "inverse_row_l2_error", "coefficient_l2_error", "selected_dimension")]
counts <- accounting[, .(registered = .N, available = sum(available)), by = experiment_id]
grouped <- coordinates[, .(squared_remainder = mean(scaled_remainder^2)), by = .(experiment_id, replicate)]
remainder <- grouped[, {
  value <- sqrt(mean(squared_remainder))
  list(metric = "Scaled influence remainder", value = value, mcse = sd(squared_remainder) / (2 * sqrt(.N) * value))
}, by = experiment_id]
curvature <- fold_replications[, .(metric = "Hessian entry error", value = mean(curvature_max_error),
  mcse = sd(curvature_max_error) / sqrt(.N)), by = experiment_id]
inverse <- fold_replications[, .(metric = "Inverse row error", value = mean(inverse_row_l2_error),
  mcse = sd(inverse_row_l2_error) / sqrt(.N)), by = experiment_id]
plot_data <- merge(rbindlist(list(remainder, curvature, inverse)),
  registry[, .(experiment_id, n_clusters, p, s)], by = "experiment_id")
plot_data[, scenario := paste0("p = ", p, "   s = ", s)]
plot_data[, `:=`(lower = value - 1.96 * mcse, upper = value + 1.96 * mcse)]
if (any(!is.finite(plot_data$value)) || any(plot_data$lower <= 0)) stop("The diagnostic logarithmic scale is unsupported.")
plot_data[, metric := factor(metric, levels = c("Hessian entry error", "Inverse row error", "Scaled influence remainder"))]
plot <- ggplot(plot_data, aes(n_clusters, value, colour = scenario, shape = scenario, linetype = scenario)) +
  geom_linerange(aes(ymin = lower, ymax = upper), linewidth = .35) +
  geom_line(linewidth = .45) + geom_point(size = 1.9) +
  facet_wrap(~metric, nrow = 1L, scales = "free_y") +
  scale_x_log10(breaks = c(100, 200, 400)) + scale_y_log10() +
  scale_colour_manual(values = c("#0072B2", "#D55E00", "#009E73", "#7B3294")) +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash")) +
  labs(x = "Independent clusters", y = "Error magnitude", colour = NULL, shape = NULL, linetype = NULL) +
  theme_bw(base_size = 9, base_family = "Times") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey90", linewidth = .2),
    legend.position = "top", legend.key.width = grid::unit(1.1, "cm"),
    strip.background = element_rect(fill = "grey96", colour = "grey60"),
    plot.margin = margin(6, 8, 6, 6)) +
  guides(colour = guide_legend(nrow = 2L), shape = guide_legend(nrow = 2L), linetype = guide_legend(nrow = 2L))
output <- "results/tables/submission_20261028/gaussian_diagnostics"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
fwrite(accounting, file.path(output, "accounting.csv"))
fwrite(coordinates, file.path(output, "coordinate_records.csv"))
fwrite(folds, file.path(output, "fold_records.csv"))
fwrite(coordinate_summary, file.path(output, "coordinate_summary.csv"))
fwrite(fold_summary, file.path(output, "fold_summary.csv"))
fwrite(plot_data, file.path(output, "figure_source.csv"))
ggsave("manuscript/generated/figures/fig_gaussian_influence_diagnostics.pdf", plot, width = 7.0, height = 3.25,
  units = "in", device = cairo_pdf)
jsonlite::write_json(list(source = root, settings = nrow(registry), registered_datasets = nrow(jobs),
  available_datasets = sum(accounting$available), statistical_unit = "independent_simulation_replication",
  uncertainty = "point_estimate_plus_or_minus_1.96_Monte_Carlo_standard_errors",
  inverse_error = "Euclidean_distance_from_the_unsmoothed_population_inverse_row",
  curvature_error = "maximum_entry_difference_at_the_full_population_target",
  within_replication_average = "requested_coordinates_estimation_folds_and_two_bandwidths",
  influence_remainder = "sqrt_n_times_paired_coefficient_error_minus_true_target_influence_average",
  figure_source = file.path(output, "figure_source.csv"), visual_inspection = FALSE,
  generated = as.character(Sys.time())), file.path(output, "manifest.json"), auto_unbox = TRUE, pretty = TRUE)
print(counts)
print(coordinate_summary)
