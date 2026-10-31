library(data.table)
library(jsonlite)

report_dir <- "results/preflight/empirical_selection_20261101/yeast_screening"
status <- read_json(file.path(report_dir, "status.json"), simplifyVector = TRUE)
stopifnot(status$complete_conditions == 44L, !status$formal_validation_used)
screening <- fread(file.path(report_dir, "comparisons.csv"))
selected <- screening[tau == 0.25 & relative_gain > 0][order(-relative_gain)]
stopifnot(nrow(selected) > 0L)
jobs <- selected[, .(trait, method = "POOLED", tau = 0.25, run_id = "tau025")]
jobs[, job_id := seq_len(.N)]
jobs[, directory := file.path("results/raw/application/yeast_20261101/discovery", trait, paste0(method, "_", run_id))]
jobs[, status := "pending"]
directory <- "results/raw/application/yeast_20261101/quantile_comparisons_20261101"
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
if (file.exists(file.path(directory, "jobs.csv"))) stop("比较清单已经存在")
fwrite(jobs, file.path(directory, "jobs.csv"))
fwrite(selected, file.path(directory, "candidate_selection.csv"))
write_json(list(stage = "discovery", selection_rule = "positive_point_gain_against_enet_at_prespecified_tau025",
                candidates = nrow(jobs), method = "exact_pooled_quantile_lasso",
                formal_validation_used = FALSE, manuscript_adoption = FALSE),
           file.path(directory, "design.json"), auto_unbox = TRUE, pretty = TRUE)
print(jobs[, .(trait, method, tau)])
