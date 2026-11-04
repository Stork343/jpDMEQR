library(data.table)
library(jsonlite)
library(yaml)

cases <- list(
  yeast = list(config = "config/discovery_yeast_20261101.yml",
    traits = fread("results/raw/application/yeast_20261101/screening_20261101/jobs.csv")[, unique(trait)]),
  yeast_genomic = list(config = "config/discovery_yeast_genomic_20261101.yml", traits = c("Ethanol", "Magnesium_Sulfate")),
  drugmatrix = list(config = "config/discovery_drugmatrix_20261101.yml", traits = "ALT"))
root <- "results/raw/application/native_mean_20261104"
for (name in names(cases)) {
  config <- read_yaml(cases[[name]]$config)
  directory <- file.path(root, name)
  if (file.exists(file.path(directory, "jobs.csv"))) stop("计算清单已经存在")
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  jobs <- data.table(trait = cases[[name]]$traits, method = "ENET", tau = NA_real_, run_id = config$mean_run_id)
  jobs[, job_id := seq_len(.N)]
  jobs[, directory := file.path(config$output_root, trait, paste0(method, "_", run_id))]
  jobs[, status := "pending"]
  fwrite(jobs, file.path(directory, "jobs.csv"))
  cat(name, nrow(jobs), "tasks\n")
}
