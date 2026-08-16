#!/usr/bin/env Rscript
full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else
  "scripts/application/04_summarise_GSE65391.R"
root <- if (file.exists("R/metrics_v2.R")) "." else
  normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", envir = environment())
source_v2_module(root, "metrics_v2", envir = environment())

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}
input_root <- get_arg("input", file.path(root, "results", "raw", "application"))
output_dir <- get_arg("output", file.path(root, "results", "tables", "application"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
metric_files <- list.files(input_root, pattern = "fold_metrics\\.csv$", recursive = TRUE, full.names = TRUE)
pred_files <- list.files(input_root, pattern = "predictions\\.csv$", recursive = TRUE, full.names = TRUE)
coef_files <- list.files(input_root, pattern = "coefficients\\.csv$", recursive = TRUE, full.names = TRUE)
inf_files <- list.files(input_root, pattern = "inference\\.csv$", recursive = TRUE, full.names = TRUE)
read_many <- function(files) if (length(files)) do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE)) else data.frame()
metrics <- read_many(metric_files); predictions <- read_many(pred_files)
coefficients <- read_many(coef_files); inference <- read_many(inf_files)
if (!nrow(metrics)) stop("No application fold_metrics.csv files found under: ", input_root)

ok <- metrics$status == "ok"
summary <- aggregate(
  cbind(pinball_loss, quantile_calibration, selected_genes, runtime_sec,
        kkt_residual, profile_identity_error, max_nuisance_gradient) ~ model + tau,
  data = metrics[ok, , drop = FALSE],
  FUN = function(x) c(mean = mean(x, na.rm = TRUE), sd = stats::sd(x, na.rm = TRUE),
                      median = stats::median(x, na.rm = TRUE), n = sum(is.finite(x)))
)
summary_flat <- data.frame(model = summary$model, tau = summary$tau, stringsAsFactors = FALSE)
for (v in setdiff(names(summary), c("model", "tau"))) {
  mat <- summary[[v]]
  for (j in seq_len(ncol(mat))) summary_flat[[paste(v, colnames(mat)[j], sep = "_")]] <- mat[, j]
}
failure <- aggregate(
  list(runs = metrics$status),
  by = list(model = metrics$model, tau = metrics$tau, status = metrics$status),
  FUN = length
)
write.csv(summary_flat, file.path(output_dir, "application_performance_summary.csv"), row.names = FALSE)
write.csv(failure, file.path(output_dir, "application_failure_summary.csv"), row.names = FALSE)
write.csv(metrics, file.path(output_dir, "application_all_fold_metrics.csv"), row.names = FALSE)
if (nrow(predictions)) write.csv(predictions, file.path(output_dir, "application_all_predictions.csv"), row.names = FALSE)
if (nrow(coefficients)) write.csv(coefficients, file.path(output_dir, "application_all_coefficients.csv"), row.names = FALSE)
if (nrow(inference)) write.csv(inference, file.path(output_dir, "application_all_inference.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("Application summaries written to: ", output_dir)
print(summary_flat)
