library(data.table)
library(jsonlite)

comparison <- fread("results/preflight/empirical_selection_20261101/yeast_quantile_comparison/decision_quantile.csv")
selected <- comparison[relative_gain > 0][order(-relative_gain)]
stopifnot(nrow(selected) > 0L, all(selected$phase == "discovery"))
jobs <- selected[, .(trait, method = "FOREST", tau = NA_real_, run_id = "tuned")]
jobs[, job_id := seq_len(.N)]
jobs[, directory := file.path("results/raw/application/yeast_20261101/discovery", trait, paste0(method, "_", run_id))]
jobs[, status := "pending"]
directory <- "results/raw/application/yeast_20261101/forest_comparisons_20261101"
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
if (file.exists(file.path(directory, "jobs.csv"))) stop("计算清单已经存在")
fwrite(jobs, file.path(directory, "jobs.csv"))
write_json(list(stage = "discovery", selection_rule = "positive_point_gain_against_enet_and_exact_quantile_lasso_at_tau025",
                formal_validation_used = FALSE, manuscript_adoption = FALSE,
                note = "candidate_selection_does_not_assert_statistical_superiority"),
           file.path(directory, "design.json"), auto_unbox = TRUE, pretty = TRUE)
print(jobs[, .(trait, method)])
