#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[[1]] else NULL
out_prefix <- if (length(args) >= 2) args[[2]] else "real_analysis"

if (!requireNamespace("jpDMEQR", quietly = TRUE)) {
  stop("Package 'jpDMEQR' is not installed. Install it first.")
}

res <- jpDMEQR::run_real_data_analysis(
  input_file = input_file,
  tau_grid = c(0.25, 0.5, 0.75),
  random_effect = "intercept_slope",
  verbose = TRUE
)

utils::write.csv(res$combined_table, paste0(out_prefix, "_combined_results.csv"), row.names = FALSE)
utils::write.csv(res$diagnostics, paste0(out_prefix, "_diagnostics.csv"), row.names = FALSE)
saveRDS(res, file = paste0(out_prefix, "_object.rds"))
message("Real-data analysis finished.")
