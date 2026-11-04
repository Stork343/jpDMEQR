library(data.table)
library(jsonlite)
library(yaml)

config <- read_yaml("config/discovery_yeast_1011_20261104.yml")
check <- read_json("results/preflight/empirical_selection_20261101/yeast_1011/input_checks.json", simplifyVector = TRUE)
stopifnot(check$status == "passed")
directory <- "results/raw/application/yeast_1011_20261104/screening_20261104"
implementation <- file.path(directory, "implementation")
if (dir.exists(directory)) stop("计算目录已经存在")
dir.create(implementation, recursive = TRUE)
files <- c("config/discovery_yeast_1011_20261104.yml", "research/empirical_selection_20261101/prediction_core.R",
           "research/yeast_1011_application_20261104/features.R", "research/yeast_1011_application_20261104/run_discovery.R")
destinations <- file.path(implementation, c("config.yml", "prediction_core.R", "features.R", "run_discovery.R"))
stopifnot(all(file.copy(files, destinations)))
jobs <- data.table(trait = "NaCl", method = c("JOINT", "ENET", "POOLED", "JOINT", "JOINT"),
                   tau = c(0.25, NA_real_, 0.25, 0.5, 0.75),
                   run_id = c("prediction_tau025", "weighted_cv", "tau025", "prediction_tau050", "prediction_tau075"))
jobs[, job_id := seq_len(.N)]
jobs[, directory := file.path(config$output_root, trait, paste0(method, "_", run_id))]
jobs[, status := "pending"]
fwrite(jobs, file.path(directory, "jobs.csv"))
print(jobs[, .(trait, method, tau)])
