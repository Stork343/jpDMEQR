library(data.table)
source <- "results/tables/submission_20261028/collections/submission_20261027_main_joint"
records <- readRDS(file.path(source, "records.rds"))
registry <- as.data.table(records$registry)
metrics <- as.data.table(records$metrics)
methods <- c("JOINT-PROFILE", "JOINT-PROFILE-SPLIT")
core <- metrics[method_id %in% methods]
groups <- fread(file.path(source, "coordinate_group_summary.csv"))[method_id %in% methods]
coordinates <- fread(file.path(source, "coordinate_summary.csv"))[method_id %in% methods]
if (nrow(registry) != 71L || nrow(core) != 45500L || nrow(groups) != 142L || nrow(coordinates) != 284L)
  stop("The complete proposed-method result accounting differs.")
if (anyDuplicated(core[, .(experiment_id, replicate)]) || anyDuplicated(core$seed))
  stop("The registered proposed-method datasets are duplicated.")
labels <- c(none = "Full predictor set", split_quantile_score = "Independent screening",
  ca_iqr_sis = "Same sample screening", variance_filter = "Variance filtering")
groups <- merge(groups, registry[, .(experiment_id, screen_method)], by = "experiment_id")
groups[, input := unname(labels[screen_method])]
if (anyNA(groups$input)) stop("An input design lacks its declared category.")
overview <- groups[, .(settings = .N, coverage = mean(coverage_all),
  coverage_mcse = sqrt(sum(coverage_mcse^2)) / .N,
  minimum_coverage = min(coverage_all), maximum_coverage = max(coverage_all),
  mean_length = mean(mean_length), length_mcse = sqrt(sum(length_mcse^2)) / .N), by = .(input, coordinate_type)]
accounting <- core[, .(registered = .N, converged = sum(optimizer_converged)), by = experiment_id]
input_counts <- merge(accounting, registry[, .(experiment_id, screen_method)], by = "experiment_id")
input_counts[, input := unname(labels[screen_method])]
input_counts <- input_counts[, .(registered = sum(registered), converged = sum(converged)), by = input]
coverage_wide <- dcast(groups, experiment_id ~ coordinate_type, value.var = c("coverage_all", "coverage_mcse", "mean_length", "length_mcse"))
coverage_wide <- merge(coverage_wide, accounting, by = "experiment_id")
overview_wide <- dcast(overview, input + settings ~ coordinate_type, value.var = c("coverage", "coverage_mcse"))
overview_wide <- merge(overview_wide, input_counts, by = "input")
overview_wide <- overview_wide[match(unname(labels), input)]
output <- "results/tables/submission_20261028/joint_study"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
fwrite(overview, file.path(output, "input_design_summary.csv"))
fwrite(accounting, file.path(output, "setting_accounting.csv"))
fwrite(core[optimizer_converged == FALSE, .(experiment_id, replicate, method_id, failure_stage, failure_message)],
  file.path(output, "numerical_failures.csv"))
generated <- "manuscript/generated/submission"
dir.create(generated, recursive = TRUE, showWarnings = FALSE)
number <- function(value) ifelse(is.na(value), "--", formatC(ifelse(abs(value) < .0005, 0, value), format = "f", digits = 3L))
uncertainty <- function(value, mcse) paste0(number(value), " (", number(mcse), ")")
write_table <- function(table, caption, label, stem, long = TRUE) {
  content <- knitr::kable(table, format = "latex", booktabs = TRUE, longtable = long, escape = FALSE,
    caption = caption, label = label, linesep = "")
  if (long) content <- kableExtra::kable_styling(content, latex_options = "repeat_header")
  writeLines(as.character(content), file.path(generated, paste0(stem, ".tex")))
}
table <- data.frame(Input = overview_wide$input, Settings = overview_wide$settings,
  Fits = paste0(overview_wide$converged, "/", overview_wide$registered),
  Active = uncertainty(overview_wide$coverage_active, overview_wide$coverage_mcse_active),
  Null = uncertainty(overview_wide$coverage_null, overview_wide$coverage_mcse_null))
write_table(table, "Full target coverage by the declared predictor input", "registered-joint-overview", "joint_overview", FALSE)
table <- data.frame(Setting = coverage_wide$experiment_id,
  Fits = paste0(coverage_wide$converged, "/", coverage_wide$registered),
  Active = uncertainty(coverage_wide$coverage_all_active, coverage_wide$coverage_mcse_active),
  Null = uncertainty(coverage_wide$coverage_all_null, coverage_wide$coverage_mcse_null))
write_table(table, "Registered full target coverage with Monte Carlo standard errors", "registered-joint-coverage", "joint_coverage")
table <- data.frame(Setting = coverage_wide$experiment_id,
  Active = uncertainty(coverage_wide$mean_length_active, coverage_wide$length_mcse_active),
  Null = uncertainty(coverage_wide$mean_length_null, coverage_wide$length_mcse_null))
write_table(table, "Registered interval lengths with Monte Carlo standard errors", "registered-joint-length", "joint_lengths")
table <- data.frame(Setting = coordinates$experiment_id,
  Coordinate = as.integer(sub("^x", "", coordinates$coordinate)),
  Available = paste0(coordinates$available_count, "/", coordinates$replications),
  Coverage = uncertainty(coordinates$coverage_all, coordinates$coverage_mcse),
  Interval = paste0("[", number(coordinates$coverage_lower), ", ", number(coordinates$coverage_upper), "]"),
  Ratio = number(coordinates$se_sd_ratio))
names(table)[6L] <- "$\\overline{s}/\\mathrm{SD}$"
write_table(table, "Coordinate coverage with exact binomial intervals and standard error ratios",
  "registered-joint-coordinates", "joint_coordinates")
table <- data.frame(Setting = registry$experiment_id, n = registry$n_clusters, p = registry$p,
  s = registry$s, tau = registry$tau, q = registry$q, Replications = registry$replications,
  Input = unname(labels[registry$screen_method]))
names(table)[2:6] <- c("$n$", "$p$", "$s$", "$\\tau$", "$q$")
write_table(table, "Registered statistical settings and dataset counts", "registered-study-design", "joint_registry")
shifts <- fread("results/tables/submission_20261028/prescreen_targets/coordinate_summary.csv")
table <- data.frame(Setting = shifts$experiment_id, Coordinate = as.integer(sub("^x", "", shifts$coordinate)),
  Shifted = paste0(shifts$shifted, "/", shifts$replications),
  Fraction = uncertainty(shifts$fraction_shifted, shifts$fraction_mcse),
  Interval = paste0("[", number(shifts$fraction_lower), ", ", number(shifts$fraction_upper), "]"),
  Displacement = uncertainty(shifts$mean_displacement, shifts$displacement_mcse))
write_table(table, "Gaussian target displacements after variance filtering and full design fitting",
  "registered-prescreen-targets", "prescreen_targets", FALSE)
get_coverage <- function(id, type) {
  row <- groups[experiment_id == id & coordinate_type == type]
  if (nrow(row) != 1L) stop("A registered coverage value is unavailable.")
  number(row$coverage_all)
}
values <- c(RegisteredJointSettings = as.character(nrow(registry)), RegisteredJointDatasets = as.character(nrow(core)),
  RegisteredJointFailures = as.character(sum(!core$optimizer_converged)),
  RegisteredJointAvailable = as.character(sum(core$optimizer_converged)),
  RegisteredFullPredictorSettings = as.character(overview[input == "Full predictor set" & coordinate_type == "active", settings]),
  RegisteredFullActiveCoverage = number(overview[input == "Full predictor set" & coordinate_type == "active", coverage]),
  RegisteredFullActiveMCSE = number(overview[input == "Full predictor set" & coordinate_type == "active", coverage_mcse]),
  RegisteredFullNullCoverage = number(overview[input == "Full predictor set" & coordinate_type == "null", coverage]),
  RegisteredFullNullMCSE = number(overview[input == "Full predictor set" & coordinate_type == "null", coverage_mcse]),
  VarianceFilterActiveCoverage = get_coverage("F05", "active"), VarianceFilterNullCoverage = get_coverage("F05", "null"),
  FullDesignActiveCoverage = get_coverage("F06", "active"), FullDesignNullCoverage = get_coverage("F06", "null"),
  VarianceFilterSecondShifted = as.character(shifts[experiment_id == "F05" & coordinate == "x00002", shifted]),
  VarianceFilterEleventhShifted = as.character(shifts[experiment_id == "F05" & coordinate == "x00011", shifted]))
writeLines(paste0("\\newcommand{\\", names(values), "}{", values, "}"), file.path(generated, "joint_values.tex"))
jsonlite::write_json(list(settings = nrow(registry), datasets = nrow(core), converged = sum(core$optimizer_converged),
  numerical_failures = sum(!core$optimizer_converged), requested_intervals = nrow(core) * 4L,
  covariance_unit = "independent_simulation_replication", input_summary_weighting = "equal_weight_over_registered_settings",
  failure_coverage = "unavailable_intervals_count_as_noncoverage", source = source, generated = as.character(Sys.time())),
  file.path(output, "manifest.json"), auto_unbox = TRUE, pretty = TRUE)
print(overview)
