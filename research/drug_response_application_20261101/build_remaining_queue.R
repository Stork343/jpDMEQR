library(data.table)
library(jsonlite)
library(yaml)

config <- read_yaml("config/discovery_gdsc_20261101.yml")
drugs <- setdiff(unlist(config$drugs), c("Dabrafenib", "HG6-64-1"))
jobs <- rbindlist(lapply(drugs, function(drug_name) data.table(
  trait = drug_name, method = c("FBIO", rep("JOINT", 3L)),
  tau = c(NA_real_, unlist(config$quantiles)),
  run_id = c("all", "qp_tau025", "qp_tau050", "qp_tau075"))))
jobs[, job_id := seq_len(.N)]
jobs[, directory := file.path("results/raw/application/drug_response_20261101/discovery", trait, paste0(method, "_", run_id))]
jobs[, status := "pending"]
directory <- "results/raw/application/drug_response_20261101/remaining_braf_20261101"
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
if (file.exists(file.path(directory, "jobs.csv"))) stop("计算清单已经存在")
fwrite(jobs, file.path(directory, "jobs.csv"))
write_json(list(stage = "discovery", drugs = drugs, source = "existing_frozen_candidate_list",
                formal_validation_used = FALSE, manuscript_adoption = FALSE),
           file.path(directory, "design.json"), auto_unbox = TRUE, pretty = TRUE)
print(jobs[, .(trait, method, tau)])
