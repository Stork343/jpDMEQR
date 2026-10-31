library(data.table)
audit <- jsonlite::read_json("results/preflight/dataset_search_20261030/milk_audit.json", simplifyVector = TRUE)
design <- jsonlite::read_json("results/preflight/dataset_search_20261030/milk_design_assessment.json", simplifyVector = TRUE)
completion <- jsonlite::read_json("results/preflight/milk_scs_20261031/completion_manifest.json", simplifyVector = TRUE)
inference <- fread("results/preflight/milk_scs_20261031/inference.csv")
if (nrow(inference) != completion$cases_registered || nrow(inference) != completion$cases_completed)
  stop("The manuscript table requires all registered cases.")
reference <- inference[setting == "reference"][order(tau)]
if (!identical(reference$tau, c(0.25, 0.5, 0.75))) stop("The reference quantiles are incomplete.")
directory <- "manuscript/generated/milk_profile"
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
format_number <- function(x, digits = 3L) sprintf(paste0("%.", digits, "f"), x)
interval <- function(lower, upper) paste0("[", format_number(lower), ",\\ ", format_number(upper), "]")
values <- c(MilkCows = audit$eligible_cows, MilkSamples = audit$eligible_observations,
  MilkRepeatCows = audit$complete_repeated_cows, MilkRepeatSamples = audit$complete_repeated_observations,
  MilkOriginalCows = audit$cows, MilkOriginalSamples = audit$observations,
  MilkInvalidCows = audit$identifier_excluded_cows, MilkInvalidSamples = audit$identifier_excluded_observations,
  MilkSpectra = audit$spectral_variables, MilkDimension = design$total_design_columns,
  MilkAOneCows = design$cow_counts_by_milk_type$cows[design$cow_counts_by_milk_type$Milk_typology == "A1"],
  MilkATwoCows = design$cow_counts_by_milk_type$cows[design$cow_counts_by_milk_type$Milk_typology == "A2"],
  MilkFits = completion$cases_completed, MilkSystems = completion$coefficient_systems,
  MilkMinSelected = completion$minimum_selected_dimension, MilkMaxSelected = completion$maximum_selected_dimension,
  MilkMedianRuntime = format_number(completion$median_runtime_seconds, 1L),
  MilkMaxKKT = sprintf("%.2e", completion$maximum_kkt),
  MilkMaxNuisance = sprintf("%.2e", completion$maximum_nuisance_gradient),
  MilkSpectralCorrelation = format_number(design$spectral_adjacent_correlation_quantiles[["50%"]], 4L))
for (j in seq_len(nrow(reference))) {
  name <- c("Lower", "Median", "Upper")[j]
  values[paste0("Milk", name, "Estimate")] <- format_number(reference$beta_tilde[j])
  values[paste0("Milk", name, "SE")] <- format_number(reference$se[j])
  values[paste0("Milk", name, "Interval")] <- interval(reference$simultaneous_lower[j], reference$simultaneous_upper[j])
}
writeLines(paste0("\\newcommand{\\", names(values), "}{", values, "}"), file.path(directory, "values.tex"))
lines <- c("\\begin{table}[t]", "\\centering", "\\caption{A2 milk profile coefficients for somatic cell score with all 1060 MIR variables available for adjustment. The pointwise intervals have 95\\% confidence level. Simultaneous intervals allocate 5\\% across the three fixed quantiles. All fits use a cow intercept and $\\Lambda=1$.}",
  "\\label{tab:milk-inference}", "\\begin{tabular}{rrrrr}", "\\toprule",
  "$\\tau$ & Estimate & Standard error & Pointwise interval & Simultaneous interval \\\\", "\\midrule")
for (i in seq_len(nrow(reference))) {
  row <- reference[i]
  lines <- c(lines, paste0(format_number(row$tau, 2L), " & ", format_number(row$beta_tilde), " & ",
    format_number(row$se), " & $", interval(row$lower, row$upper), "$ & $",
    interval(row$simultaneous_lower, row$simultaneous_upper), "$ \\\\"))
}
writeLines(c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}"), file.path(directory, "reference.tex"))
labels <- c(reference = "Reference", paired = "Two milk controls", bandwidth_075 = "$0.75h$",
  bandwidth_125 = "$1.25h$", penalty_050 = "$\\Lambda=0.5$", penalty_200 = "$\\Lambda=2$")
lines <- c("\\begin{table}[p]", "\\centering", "\\small",
  "\\caption{Prespecified sensitivity analyses. Each setting reports its own population profile parameter. Intervals allocate 5\\% across its three quantiles. Selected dimensions are the minimum and maximum over the two estimation samples. Time is the complete coordinate-fitting time in seconds after the shared design regression has been computed.}",
  "\\label{tab:milk-sensitivity}", "\\begin{tabular}{lrrrrr}", "\\toprule",
  "Setting & $\\tau$ & Estimate & Simultaneous interval & Dimension & Seconds \\\\", "\\midrule")
for (i in seq_len(nrow(inference))) {
  row <- inference[i]
  lines <- c(lines, paste0(labels[[row$setting]], " & ", format_number(row$tau, 2L), " & ",
    format_number(row$beta_tilde), " & $", interval(row$simultaneous_lower, row$simultaneous_upper), "$ & ",
    row$dimension_min, "--", row$dimension_max, " & ", format_number(row$runtime_seconds, 1L), " \\\\"))
}
writeLines(c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}"), file.path(directory, "sensitivity.tex"))
cat("Generated milk application values and tables from complete fitted results.\n")
