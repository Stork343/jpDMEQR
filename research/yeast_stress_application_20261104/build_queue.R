library(data.table)
library(jsonlite)
library(yaml)

config <- read_yaml("config/discovery_yeast_stress_20261104.yml")
check <- read_json("results/preflight/empirical_selection_20261101/yeast_stress/input_checks.json", simplifyVector = TRUE)
stopifnot(check$status == "passed")
directory <- "results/raw/application/yeast_stress_20261104/screening_20261104"
implementation <- file.path(directory, "implementation")
if (dir.exists(directory)) stop("计算目录已经存在")
dir.create(implementation, recursive = TRUE)
files <- c("config/discovery_yeast_stress_20261104.yml", "research/empirical_selection_20261101/prediction_core.R",
           "research/yeast_1011_application_20261104/features.R", "research/yeast_stress_application_20261104/features.R",
           "research/yeast_stress_application_20261104/run_discovery.R")
destinations <- file.path(implementation, c("config.yml", "prediction_core.R", "genomic_features.R", "features.R", "run_discovery.R"))
stopifnot(all(file.copy(files, destinations)))
jobs <- rbindlist(lapply(names(config$assays), function(assay_name)
  data.table(trait = assay_name, method = c("JOINT", "ENET"), tau = c(0.25, NA_real_),
             run_id = c("prediction_tau025", "weighted_cv"))))
jobs[, job_id := seq_len(.N)]
jobs[, directory := file.path(config$output_root, trait, paste0(method, "_", run_id))]
jobs[, status := "pending"]
fwrite(jobs, file.path(directory, "jobs.csv"))
print(jobs[, .(trait, method, tau)])
