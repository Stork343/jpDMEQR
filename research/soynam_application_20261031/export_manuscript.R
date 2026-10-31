library(data.table)
library(yaml)
library(jsonlite)
config <- read_yaml("config/application_soynam_20261031.yml")
arguments <- commandArgs(trailingOnly = TRUE)
mode <- if (length(arguments)) arguments[1L] else "all"
if (!mode %in% c("all", "joint")) stop("Unknown manuscript export mode.")
source <- "results/tables/soynam_20261031"
prediction_source <- if (mode == "all") source else "results/preflight/soynam_20261031/joint_prediction"
destination <- "manuscript/generated/soynam_profile"
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
effects <- fread(file.path(source, "inference.csv"))
contrasts <- fread(file.path(source, "reference_contrasts.csv"))
all_contrasts <- fread(file.path(source, "all_setting_contrasts.csv"))
diagnostics <- fread(file.path(source, "diagnostics.csv"), integer64 = "double")
metrics <- fread(file.path(prediction_source, "prediction_metrics.csv"))
comparisons <- fread(file.path(prediction_source, "prediction_comparisons.csv"))
year_metrics <- fread(file.path(prediction_source, "prediction_year_metrics.csv"))
descriptive <- fread(file.path(source, "descriptive.csv"))
audit <- read_json("results/preflight/soynam_20261031/data_audit.json", simplifyVector = TRUE)
manifest <- read_json("results/preflight/soynam_20261031/input_manifest.json", simplifyVector = TRUE)
input <- readRDS(file.path(config$results_directory, "input.rds"))
expected_methods <- if (mode == "all") 3L else 2L
if (nrow(effects) != 42L || nrow(contrasts) != 11L || nrow(metrics) != 3L * expected_methods ||
    nrow(comparisons) != 6L * (expected_methods - 1L))
  stop("The complete empirical table set is required for manuscript export.")
fmt <- function(x, digits = 2L) sprintf(paste0("%.", digits, "f"), x)
interval <- function(lower, upper, digits = 2L) paste0("$[", fmt(lower, digits), ",\\ ", fmt(upper, digits), "]$")
parameter <- function(name) {
  found <- contrasts[parameter == name]
  if (nrow(found) != 1L) stop("A manuscript parameter is not unique.")
  found
}
values <- list(Lines = manifest$lines, Families = audit$families, Records = manifest$records,
  OriginalMarkers = audit$markers, OriginalLines = audit$lines, OriginalRecords = audit$observations,
  Markers = manifest$genomic_markers, Dimension = manifest$total_columns, CheckPlots = manifest$check_records,
  TestLines = sum(input$line_folds$prediction_fold == config$prediction$heldout_fold),
  TrainingLines = sum(input$line_folds$prediction_fold != config$prediction$heldout_fold),
  YieldMeanThirteen = fmt(descriptive[Year == 2013 & variable == "YLD", mean]),
  YieldMeanFourteen = fmt(descriptive[Year == 2014 & variable == "YLD", mean]),
  Fits = uniqueN(effects$case_id), MinSelected = min(diagnostics$dimension),
  MaxSelected = max(diagnostics$dimension), MaxKKT = sprintf("%.2e", max(diagnostics$kkt)),
  MaxNuisance = sprintf("%.2e", max(diagnostics$nuisance_gradient)),
  ReferenceMinMinutes = fmt(min(effects[setting == "reference", runtime_seconds]) / 60, 1),
  ReferenceMaxMinutes = fmt(max(effects[setting == "reference", runtime_seconds]) / 60, 1))
mapping <- c(MedianThirteen = "ACC10_2013_q50", MedianFourteen = "ACC10_2014_q50",
  UpperYearDifference = "year_2014_minus_2013_q75", LowerYearDifference = "year_2014_minus_2013_q25",
  MedianYearDifference = "year_2014_minus_2013_q50")
for (name in names(mapping)) {
  row <- parameter(mapping[[name]])
  values[[name]] <- fmt(row$estimate)
  values[[paste0(name, "Interval")]] <- interval(row$lower, row$upper)
}
reference <- effects[setting == "reference", .(coordinate, tau, reference = estimate_bu_ac)]
joined <- merge(effects, reference, by = c("coordinate", "tau"))
values$BandwidthChange <- fmt(joined[setting %in% c("bandwidth_075", "bandwidth_125"),
  max(abs(estimate_bu_ac - reference))], 3)
values$CoefficientMinimum <- fmt(min(reference$reference))
values$CoefficientMaximum <- fmt(max(reference$reference))
year_contrasts <- all_contrasts[parameter == "year_2014_minus_2013_q75"]
if (nrow(year_contrasts) != 7L) stop("The sensitivity contrast family is incomplete.")
values$YearContrastMinimum <- fmt(min(year_contrasts$estimate))
values$YearContrastMaximum <- fmt(max(year_contrasts$estimate))
spatial_contrast <- year_contrasts[setting == "spatial_polynomial"]
values$SpatialYearContrast <- fmt(spatial_contrast$estimate)
values$SpatialYearContrastInterval <- interval(spatial_contrast$lower, spatial_contrast$upper)
labels <- c(joint_canopy = "Joint profile with canopy", joint_genomic = "Joint profile without canopy",
  pooled_quantile_lasso = "Pooled quantile Lasso")
quantile_names <- c("Lower", "Median", "Upper")
for (q in seq_along(config$quantiles)) for (comparator_id in unique(comparisons$comparator)) {
  row <- comparisons[tau == config$quantiles[[q]] & comparator == comparator_id & calibrated == TRUE]
  if (nrow(row) != 1L) stop("A prediction comparison is not unique.")
  name <- paste0(if (comparator_id == "joint_genomic") "Genomic" else "Pooled", quantile_names[q])
  values[[paste0(name, "Gain")]] <- fmt(row$relative_percent, 1)
  values[[paste0(name, "GainInterval")]] <- interval(row$relative_lower, row$relative_upper, 1)
  values[[paste0(name, "LossDifference")]] <- fmt(row$improvement, 3)
}
writeLines(vapply(names(values), function(name) paste0("\\newcommand{\\Soy", name, "}{", values[[name]], "}"), character(1)),
  file.path(destination, "values.tex"))
table <- function(filename, caption, label, alignment, header, rows, size = "\\small") {
  writeLines(c("\\begin{table}[t]", "\\centering", size, paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"), paste0("\\begin{tabular}{", alignment, "}"), "\\toprule",
    paste0(header, " \\\\"), "\\midrule", paste0(rows, " \\\\"), "\\bottomrule", "\\end{tabular}", "\\end{table}"),
    file.path(destination, paste0(filename, ".tex")))
}
ref <- effects[setting == "reference"]
setorder(ref, coordinate, tau)
rows <- vapply(seq_len(nrow(ref)), function(i) paste(sub("ACC10_", "", ref$coordinate[i]), fmt(ref$tau[i]),
  fmt(ref$estimate_bu_ac[i], 3), fmt(ref$se_bu_ac[i], 3),
  interval(ref$simultaneous_lower[i], ref$simultaneous_upper[i], 3), sep = " & "), character(1))
table("reference", "Canopy profile coefficients in bu/ac per ten percentage points of coverage. Simultaneous intervals allocate 5\\% across six coordinates and five contrasts. The reference penalty is $\\Lambda=1$.",
  "tab:soy-reference", "rrrrr", "Year & $\\tau$ & Estimate & Standard error & Simultaneous interval", rows)
contrast_names <- c(upper_minus_lower_2013 = "Upper minus lower in 2013",
  upper_minus_lower_2014 = "Upper minus lower in 2014",
  year_2014_minus_2013_q25 = "2014 minus 2013 at $\\tau=0.25$",
  year_2014_minus_2013_q50 = "2014 minus 2013 at $\\tau=0.50$",
  year_2014_minus_2013_q75 = "2014 minus 2013 at $\\tau=0.75$")
contrast <- contrasts[parameter %in% names(contrast_names)]
rows <- vapply(seq_len(nrow(contrast)), function(i) paste(contrast_names[[contrast$parameter[i]]],
  fmt(contrast$estimate[i], 3), fmt(contrast$se[i], 3), interval(contrast$lower[i], contrast$upper[i], 3), sep = " & "), character(1))
table("contrasts", "Prespecified canopy coefficient differences in bu/ac per ten percentage points. These intervals use the same eleven-parameter family as the coordinate intervals.",
  "tab:soy-contrasts", "lrrr", "Contrast & Estimate & Standard error & Simultaneous interval", rows)
descriptive_names <- c(YLD = "Yield (bu/ac)", ACC = "Canopy coverage (\\%)", DTM = "Maturity (days)")
rows <- vapply(seq_len(nrow(descriptive)), function(i) paste(descriptive$Year[i],
  descriptive_names[[descriptive$variable[i]]], fmt(descriptive$mean[i]), fmt(descriptive$sd[i]),
  fmt(descriptive$q25[i]), fmt(descriptive$median[i]), fmt(descriptive$q75[i]), sep = " & "), character(1))
table("descriptive", paste0("Observed distributions in the complete two-year panel. Each year contains one observation for each of the same ", manifest$lines, " lines."),
  "tab:soy-descriptive", "rlrrrrr", "Year & Variable & Mean & SD & Quartile 1 & Median & Quartile 3", rows)
setting_names <- c(reference = "Reference", bandwidth_075 = "$0.75h$", bandwidth_125 = "$1.25h$",
  penalty_050 = "$\\Lambda=0.5$", penalty_200 = "$\\Lambda=2$",
  maturity_records = "Maturity restriction", spatial_polynomial = "Spatial polynomial")
rows <- character()
for (setting_id in names(setting_names)) for (tau_value in unlist(config$quantiles)) {
  current <- effects[setting == setting_id & tau == tau_value]
  if (nrow(current) != 2L) stop("A sensitivity setting lacks its two coordinates.")
  setorder(current, coordinate)
  rows <- c(rows, paste(setting_names[[setting_id]], fmt(tau_value),
    paste0(fmt(current$estimate_bu_ac[1L], 3), " ", interval(current$simultaneous_lower[1L], current$simultaneous_upper[1L], 3)),
    paste0(fmt(current$estimate_bu_ac[2L], 3), " ", interval(current$simultaneous_lower[2L], current$simultaneous_upper[2L], 3)), sep = " & "))
}
table("sensitivity", "Complete registered sensitivity calculations. Entries are estimates and simultaneous intervals in bu/ac per ten percentage points. Each setting uses the eleven-parameter critical value; the working penalty belongs to its profile target.",
  "tab:soy-sensitivity", "lrrr", "Setting & $\\tau$ & 2013 & 2014", rows, "\\footnotesize")
rows <- vapply(seq_len(nrow(year_contrasts)), function(i) paste(setting_names[[year_contrasts$setting[i]]],
  fmt(year_contrasts$estimate[i], 3), fmt(year_contrasts$se[i], 3),
  interval(year_contrasts$lower[i], year_contrasts$upper[i], 3), sep = " & "), character(1))
table("contrast_sensitivity", "The 2014 minus 2013 canopy coefficient contrast at the upper quartile across all registered settings. Units are bu/ac per ten percentage points of coverage; intervals use the eleven-parameter critical value within each setting.",
  "tab:soy-contrast-sensitivity", "lrrr", "Setting & Estimate & Standard error & Simultaneous interval", rows)
rows <- vapply(unlist(config$quantiles), function(tau_value) {
  current <- metrics[tau == tau_value]
  cells <- c(fmt(tau_value), fmt(current[method == "joint_canopy", loss], 3),
    fmt(current[method == "joint_genomic", loss], 3))
  if (mode == "all") cells <- c(cells, fmt(current[method == "pooled_quantile_lasso", loss], 3))
  paste(cells, collapse = " & ")
}, character(1))
table("prediction", "Mean pinball loss on the same held-out lines, in bu/ac. Each line receives equal weight across its two observations. Every predictor uses calibration estimated from training validation residuals.",
  "tab:soy-prediction", if (mode == "all") "rrrr" else "rrr",
  paste0("$\\tau$ & Joint with canopy & Joint without canopy", if (mode == "all") " & Pooled with canopy" else ""), rows)
if (mode == "all") {
  pooled <- comparisons[comparator == "pooled_quantile_lasso" & calibrated == TRUE]
  rows <- vapply(seq_len(nrow(pooled)), function(i) paste(fmt(pooled$tau[i]),
    fmt(pooled$improvement[i], 3), interval(pooled$lower[i], pooled$upper[i], 3),
    fmt(pooled$relative_percent[i], 2), interval(pooled$relative_lower[i], pooled$relative_upper[i], 2), sep = " & "), character(1))
  table("pooled_comparison", "Held-out prediction comparison between joint profile and pooled quantile Lasso, both with canopy and genomic variables. Positive improvements favour the joint predictor. Absolute values have bu/ac units, and relative values are percentages. Intervals allocate 5\\% across the six registered method and quantile comparisons.",
    "tab:soy-pooled-comparison", "rrrrr", "$\\tau$ & Improvement & Simultaneous interval & Relative improvement & Simultaneous interval", rows)
}
rows <- vapply(seq_len(nrow(metrics)), function(i) paste(labels[[metrics$method[i]]], fmt(metrics$tau[i]),
  fmt(metrics$raw_loss[i], 3), fmt(metrics$loss[i], 3), fmt(metrics$raw_calibration[i], 3),
  fmt(metrics$calibration[i], 3), sep = " & "), character(1))
table("prediction_details", "Raw and training-calibrated predictions on held-out lines. Coverage is the fraction of observed yields at or below the predicted quantile. The nominal fractions are the corresponding values of $\\tau$.",
  "tab:soy-prediction-details", "lrrrrr", "Method & $\\tau$ & Raw loss & Calibrated loss & Raw coverage & Calibrated coverage", rows, "\\footnotesize")
writeLines(c("\\begin{figure}[t]", "\\centering",
  "\\includegraphics[width=\\textwidth]{generated/soynam_profile/results.pdf}",
  "\\caption{Canopy information in the soybean trials. Panel A gives the year-specific profile coefficients per ten percentage points of coverage, with 95\\% Bonferroni simultaneous intervals for the eleven prespecified coordinate and contrast parameters. Panel B gives the percentage reduction in held-out pinball loss when canopy is added to the joint predictor using genomic and field information. Its intervals use paired complete-line losses and the registered six-comparison Bonferroni adjustment. The prediction comparison uses training-only intercept calibration.}",
  "\\label{fig:soy-results}", "\\end{figure}"), file.path(destination, "figure.tex"))
cat("Generated manuscript values and tables for mode:", mode, "\n")
