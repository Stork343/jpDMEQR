library(data.table)
library(jsonlite)
library(yaml)

config <- read_yaml("config/development_ridge_offset_forest_20261105.yml")
screen <- fread(config$screening_decision)
stopifnot(setequal(screen[further_comparisons_required == TRUE]$assay, unlist(config$assays)))
check <- read_json("results/preflight/empirical_selection_20261101/ridge_direction/forest_checks.json")
stopifnot(check$status == "passed")
root <- config$output_root
if (dir.exists(root)) stop("比较目录已经存在")
implementation <- file.path(root, "implementation")
dir.create(implementation, recursive = TRUE)
files <- c("config/development_ridge_offset_forest_20261105.yml",
  "research/empirical_selection_20261101/prediction_core.R", "research/empirical_selection_20261101/ridge_direction.R",
  "research/yeast_stress_application_20261104/features.R", "research/yeast_proteome_application_20261105/features.R",
  "research/ridge_direction_20261105/forest.R", "research/ridge_direction_20261105/run_forest.R")
destinations <- file.path(implementation, c("config.yml", "prediction_core.R", "ridge_direction.R", "contrasts.R",
                                           "features.R", "forest.R", "run_discovery.R"))
stopifnot(all(file.copy(files, destinations)))
jobs <- rbindlist(lapply(unlist(config$assays), function(assay_name)
  data.table(trait = assay_name, method = c("QRF_RAW", "QRF_OFFSET"), tau = NA_real_, run_id = "native")))
jobs[, job_id := seq_len(.N)]
jobs[, directory := file.path(root, trait, paste0(method, "_", run_id))]
jobs[, status := "pending"]
fwrite(jobs, file.path(root, "jobs.csv"))
print(jobs[, .(trait, method)])
