library(data.table)
library(jsonlite)
library(yaml)

config <- read_yaml("config/discovery_drugmatrix_20261101.yml")
check <- read_json("results/preflight/empirical_selection_20261101/drugmatrix_metadata/input_checks.json", simplifyVector = TRUE)
stopifnot(check$status == "passed")
jobs <- data.table(trait = "ALT", method = c("JOINT", "ENET", "POOLED", "JOINT", "JOINT"),
                   tau = c(0.75, NA_real_, 0.75, 0.25, 0.5),
                   run_id = c("qp_tau075", "cv", "tau075", "qp_tau025", "qp_tau050"))
jobs[, job_id := seq_len(.N)]
jobs[, directory := file.path(config$output_root, trait, paste0(method, "_", run_id))]
jobs[, status := "pending"]
directory <- "results/raw/application/drugmatrix_20261101/screening_20261101"
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
if (file.exists(file.path(directory, "jobs.csv"))) stop("计算清单已经存在")
fwrite(jobs, file.path(directory, "jobs.csv"))
print(jobs[, .(trait, method, tau)])
