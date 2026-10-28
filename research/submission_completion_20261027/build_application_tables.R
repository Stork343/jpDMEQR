library(data.table)
source <- "results/tables/submission_20261028/application"
manifest <- jsonlite::read_json(file.path(source, "manifest.json"), simplifyVector = TRUE)
if (manifest$registered_comparator_jobs != 1350L || manifest$registered_joint_jobs != 900L)
  stop("The completed application study is required.")
summary <- fread(file.path(source, "paired_comparison_summary.csv"))
accounting <- fread(file.path(source, "job_accounting.csv"))
output <- "manuscript/generated/submission"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
number <- function(x) ifelse(is.na(x), "--", formatC(ifelse(abs(x) < .0005, 0, x), format = "f", digits = 3L))
labels <- c("POOLED-QR-LASSO" = "Pooled Lasso", "QGEE-SCAD" = "QGEE SCAD", LQMM = "LQMM")
for (processing in unique(summary$processing)) {
  for (scope in c("available_repetitions", "five_repetitions")) {
    selected <- summary[summary$processing == processing & summary$scope == scope]
    setorder(selected, model, tau, method)
    if (nrow(selected) != 18L) stop("An application table lacks a method setting.")
    table <- data.frame(Model = selected$model, Quantile = number(selected$tau),
      Method = unname(labels[selected$method]), Fits = paste0(selected$successful_fits, "/", selected$attempted_fits),
      Subjects = selected$subjects,
      Repetitions = ifelse(selected$minimum_repetitions == selected$maximum_repetitions,
        selected$minimum_repetitions, paste0(selected$minimum_repetitions, "--", selected$maximum_repetitions)),
      Difference = number(selected$difference),
      Range = ifelse(is.na(selected$lower), "--", paste0("[", number(selected$lower), ", ", number(selected$upper), "]")))
    names(table) <- c("Model", "$\\tau$", "Method", "Fits", "Subjects", "Repeats", "$\\Delta L$", "Range")
    route <- switch(processing, reference = "reference processing",
      official_R1_R2_normexp_training_quantile = "normalised raw arrays",
      processed_median_across_mapped_probes = "median probe aggregation", stop("Unknown processing route."))
    suffix <- if (scope == "available_repetitions") "available" else "complete"
    caption <- paste0("Paired subject prediction comparisons under ", route,
      if (scope == "available_repetitions") " using available repetitions" else " using subjects with all five repetitions")
    stem <- paste0("application_", if (processing == "reference") "reference" else if (grepl("official", processing)) "raw" else "median", "_", suffix)
    content <- knitr::kable(table, format = "latex", booktabs = TRUE, longtable = TRUE,
      caption = caption, label = stem, escape = FALSE, linesep = "", align = c("l", "r", "l", rep("r", 5L)))
    writeLines(as.character(content), file.path(output, paste0(stem, ".tex")))
    if (scope == "available_repetitions") {
      performance <- data.frame(Model = selected$model, Quantile = number(selected$tau),
        Method = unname(labels[selected$method]), Fits = paste0(selected$successful_fits, "/", selected$attempted_fits),
        Subjects = selected$subjects, Loss = number(selected$loss), Calibration = number(selected$calibration))
      names(performance)[2L] <- "$\\tau$"
      performance_stem <- sub("_available$", "_performance", stem)
      content <- knitr::kable(performance, format = "latex", booktabs = TRUE, longtable = TRUE,
        caption = paste0("Available subject prediction losses and calibration under ", route),
        label = performance_stem, escape = FALSE, linesep = "", align = c("l", "r", "l", rep("r", 4L)))
      writeLines(as.character(content), file.path(output, paste0(performance_stem, ".tex")))
    }
  }
}
counts <- accounting[, .(attempted = .N, successful = sum(successful), failed = sum(!successful)), by = method]
fwrite(counts, file.path(source, "method_accounting.csv"))
macros <- c(paste0("\\newcommand{\\SubmissionComparatorFits}{", nrow(accounting), "}"),
  paste0("\\newcommand{\\SubmissionComparatorSuccessful}{", sum(accounting$successful), "}"),
  paste0("\\newcommand{\\SubmissionComparatorFailed}{", sum(!accounting$successful), "}"))
writeLines(macros, file.path(output, "application_values.tex"))
cat("Generated nine complete application comparison tables\n")
