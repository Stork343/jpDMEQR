library(data.table)
library(jsonlite)

directory <- "results/raw/application/yeast_20261101/quantile_comparisons_20261101"
jobs <- fread(file.path(directory, "progress.csv"))
stopifnot(all(jobs$status == "completed"))
tables <- list()
for (trait_name in jobs$trait) {
  log_path <- file.path(directory, "logs", paste0("report_", trait_name, ".log"))
  code <- system2(file.path(R.home("bin"), "Rscript"),
    c("research/yeast_application_20261101/report_comparators.R",
      "config/discovery_yeast_20261101.yml", shQuote(trait_name)),
    stdout = log_path, stderr = log_path)
  stopifnot(code == 0L)
  tables[[trait_name]] <- fread(file.path("results/preflight/empirical_selection_20261101/yeast_comparators",
    "yeast_20261101", trait_name, "comparisons.csv"))
}
table <- rbindlist(tables)
output <- "results/preflight/empirical_selection_20261101/yeast_quantile_comparison"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
fwrite(table, file.path(output, "comparisons.csv"))
decision <- table[tau == 0.25 & comparator == "POOLED"]
fwrite(decision[order(-relative_gain)], file.path(output, "decision_quantile.csv"))
write_json(list(stage = "discovery", conditions = nrow(decision),
                positive_point_gain = sum(decision$relative_gain > 0),
                pointwise_lower_bound_above_zero = sum(decision$lower > 0),
                formal_validation_used = FALSE, manuscript_adoption = FALSE),
           file.path(output, "status.json"), auto_unbox = TRUE, pretty = TRUE)
print(decision[order(-relative_gain)])
