library(data.table)
library(jsonlite)
library(yaml)

config <- read_yaml("config/discovery_arabidopsis_20261105.yml")
check <- read_json("results/preflight/empirical_selection_20261101/arabidopsis/input_checks.json", simplifyVector = TRUE)
stopifnot(check$status == "passed")
solver_check <- read_json("results/preflight/empirical_selection_20261101/arabidopsis/solver_checks.json", simplifyVector = TRUE)
stopifnot(solver_check$status == "passed")
directory <- "results/raw/application/arabidopsis_20261105/screening_20261105"
implementation <- file.path(directory, "implementation")
if (dir.exists(directory)) stop("计算目录已经存在")
dir.create(implementation, recursive = TRUE)
files <- c("config/discovery_arabidopsis_20261105.yml", "research/empirical_selection_20261101/prediction_core.R",
           "research/arabidopsis_application_20261105/features.R", "research/arabidopsis_application_20261105/run_discovery.R")
destinations <- file.path(implementation, c("config.yml", "prediction_core.R", "features.R", "run_discovery.R"))
stopifnot(all(file.copy(files, destinations)))
jobs <- data.table(trait = "Flowering", method = c("JOINT", "ENET", "POOLED", "JOINT", "JOINT"),
                   tau = c(0.25, NA_real_, 0.25, 0.5, 0.75),
                   run_id = c("prediction_tau025", "weighted_cv", "tau025", "prediction_tau050", "prediction_tau075"))
jobs[, job_id := seq_len(.N)]
jobs[, directory := file.path(config$output_root, trait, paste0(method, "_", run_id))]
jobs[, status := "pending"]
fwrite(jobs, file.path(directory, "jobs.csv"))
print(jobs[, .(trait, method, tau)])
