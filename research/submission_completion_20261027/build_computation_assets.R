library(data.table)
root <- "results/raw/simulation/submission_20261028_computation"
if (!file.exists(file.path(root, "completed.rds"))) stop("The registered computation study is incomplete.")
source <- readRDS("results/tables/submission_20261028/collections/submission_20261028_computation/records.rds")
metrics <- as.data.table(source$metrics)
registry <- as.data.table(source$registry)
request <- readRDS(file.path(root, "request.rds"))
if (nrow(request$jobs) != 350L || nrow(metrics) != 1000L || any(!metrics$optimizer_converged))
  stop("The complete computation accounting differs.")
parallel_checks <- list()
for (replicate in 1:50) {
  first <- readRDS(file.path(root, "replications/G06", sprintf("replicate_%05d/result.rds", replicate)))
  second <- readRDS(file.path(root, "replications/G07", sprintf("replicate_%05d/result.rds", replicate)))
  if (first$seed != second$seed || !identical(first$computational_seed_id, second$computational_seed_id))
    stop("Parallel computation seeds differ.")
  for (method in first$metrics$method_id) {
    first_coordinate <- first$coordinates[first$coordinates$method_id == method, ]
    second_coordinate <- second$coordinates[second$coordinates$method_id == method, ]
    if (!identical(first_coordinate$coordinate, second_coordinate$coordinate)) stop("Parallel target coordinates differ.")
    beta_first <- first$method_diagnostics[[method]]$beta_hat
    beta_second <- second$method_diagnostics[[method]]$beta_hat
    if (!identical(names(beta_first), names(beta_second))) stop("Parallel prediction coordinates differ.")
    error <- max(abs(beta_first - beta_second), 0)
    for (field in c("beta_tilde", "estimated_se", "ci_lower", "ci_upper")) {
      left <- first_coordinate[[field]]
      right <- second_coordinate[[field]]
      if (!identical(is.na(left), is.na(right))) stop("Parallel interval availability differs.")
      error <- max(error, abs(left[!is.na(left)] - right[!is.na(right)]), 0)
    }
    if (!is.finite(error) || error > 1e-8) stop("Parallel fitted values differ beyond numerical tolerance.")
    parallel_checks[[length(parallel_checks) + 1L]] <- data.table(replicate = replicate, method_id = method,
      maximum_difference = error)
  }
}
timings <- metrics[, .(replications = .N, converged = sum(optimizer_converged),
  median_seconds = median(runtime_sec), lower_quartile_seconds = unname(quantile(runtime_sec, .25)),
  upper_quartile_seconds = unname(quantile(runtime_sec, .75)), maximum_seconds = max(runtime_sec)),
  by = .(experiment_id, method_id)]
parts <- vector("list", nrow(request$jobs))
for (i in seq_len(nrow(request$jobs))) {
  job <- request$jobs[i, ]
  result <- readRDS(file.path(root, "replications", job$experiment_id, sprintf("replicate_%05d/result.rds", job$replicate)))
  timing <- result$method_diagnostics[["JOINT-PROFILE"]]$submission$fit$timing
  if (!all(c("training", "inference", "total") %in% names(timing)) || any(!is.finite(timing)))
    stop("A proposed-method component clock is unavailable.")
  parts[[i]] <- data.table(experiment_id = job$experiment_id, replicate = job$replicate,
    selection_seconds = timing[["training"]], inference_seconds = timing[["inference"]], total_seconds = timing[["total"]])
}
parts <- rbindlist(parts)
part_summary <- parts[, lapply(.SD, median), by = experiment_id,
  .SDcols = c("selection_seconds", "inference_seconds", "total_seconds")]
wall <- fread(file.path(root, "batch_wall_times.csv"))
output <- "results/tables/submission_20261028/computation"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
fwrite(timings, file.path(output, "timing_summary.csv"))
fwrite(parts, file.path(output, "component_records.csv"))
fwrite(part_summary, file.path(output, "component_summary.csv"))
fwrite(rbindlist(parallel_checks), file.path(output, "parallel_reproduction.csv"))
fwrite(wall, file.path(output, "batch_wall_times.csv"))
generated <- "manuscript/generated/submission"
dir.create(generated, recursive = TRUE, showWarnings = FALSE)
number <- function(x) formatC(x, format = "f", digits = 3L)
labels <- c("JOINT-PROFILE" = "Joint profile", "JOINT-PROFILE-TRUE-SUPPORT" = "Known support oracle", "POOLED-QR-LASSO" = "Pooled Lasso")
available <- metrics[method_id == "JOINT-PROFILE", .(lower = min(p), upper = max(p)), by = experiment_id]
available <- available[match(registry$experiment_id, available$experiment_id)]
design <- data.frame(Setting = registry$experiment_id, n = registry$n_clusters, p = registry$p,
  Available = ifelse(available$lower == available$upper, as.character(available$lower),
    paste0(available$lower, "--", available$upper)),
  q = registry$q, Targets = registry$target_coordinate_count, Workers = registry$workers)
names(design)[c(2L, 3L, 5L)] <- c("$n$", "$p$", "$q$")
writeLines(as.character(knitr::kable(design, format = "latex", booktabs = TRUE, escape = FALSE,
  caption = "Registered computation settings with fifty timed repetitions", label = "registered-computation-design", linesep = "")),
  file.path(generated, "computation_design.tex"))
table <- data.frame(Setting = timings$experiment_id, Method = unname(labels[timings$method_id]),
  Fits = paste0(timings$converged, "/", timings$replications), Median = number(timings$median_seconds),
  Quartiles = paste0("[", number(timings$lower_quartile_seconds), ", ", number(timings$upper_quartile_seconds), "]"))
writeLines(as.character(knitr::kable(table, format = "latex", booktabs = TRUE, longtable = TRUE,
  caption = "Elapsed fitting time in seconds including tuning", label = "registered-computation-time", linesep = "")),
  file.path(generated, "computation_time.tex"))
values <- c(SubmissionTimingReplications = as.character(nrow(request$jobs)),
  SubmissionTimingJointHigh = number(timings[experiment_id == "G02" & method_id == "JOINT-PROFILE", median_seconds]),
  SubmissionTimingPooledHigh = number(timings[experiment_id == "G02" & method_id == "POOLED-QR-LASSO", median_seconds]),
  SubmissionTimingTwelve = number(timings[experiment_id == "G05" & method_id == "JOINT-PROFILE", median_seconds]),
  SubmissionWallFour = number(wall[experiment_id == "G06", wall_seconds]),
  SubmissionWallEight = number(wall[experiment_id == "G07", wall_seconds]),
  SubmissionParallelComparisons = as.character(length(parallel_checks)))
writeLines(paste0("\\newcommand{\\", names(values), "}{", values, "}"), file.path(generated, "computation_values.tex"))
jsonlite::write_json(list(registered_timed_runs = nrow(request$jobs), method_fits = nrow(metrics),
  replication_unit = "one_generated_dataset_and_requested_method_set", method_clock = "tuning_and_fitting_after_prescreening_excluding_data_generation",
  throughput_clock = "complete_batch_including_worker_startup_and_data_generation", timing_summary = "median_and_interquartile_range",
  parallel_parameter_comparisons = length(parallel_checks),
  maximum_parallel_difference = max(vapply(parallel_checks, function(x) x$maximum_difference, numeric(1))),
  source = root, implementation_commit = request$implementation_commit,
  generated = as.character(Sys.time())), file.path(output, "manifest.json"), auto_unbox = TRUE, pretty = TRUE)
print(timings)
print(wall)
