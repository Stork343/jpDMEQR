library(data.table)
library(yaml)

cases <- list(rice = list(config = "config/discovery_rice_baap_20261104.yml", trait = "Cadmium"),
               drugmatrix = list(config = "config/discovery_drugmatrix_20261101.yml", trait = "ALT"))
for (name in names(cases)) {
  config <- read_yaml(cases[[name]]$config)
  jobs <- data.table(trait = cases[[name]]$trait, method = c("JOINT", "ENET", "JOINT", "JOINT"),
                     tau = c(0.75, NA_real_, 0.25, 0.5))
  jobs[, run_id := ifelse(method == "JOINT",
    paste0(config$joint_run_prefix, "_tau", sprintf("%03d", as.integer(round(tau * 100)))), config$mean_run_id)]
  jobs[, job_id := seq_len(.N)]
  jobs[, directory := file.path(config$output_root, trait, paste0(method, "_", run_id))]
  jobs[, status := "pending"]
  directory <- file.path("results/raw/application/prediction_20261104", name)
  if (file.exists(file.path(directory, "jobs.csv"))) stop("计算清单已经存在")
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  fwrite(jobs, file.path(directory, "jobs.csv"))
  print(jobs[, .(trait, method, tau, run_id)])
}
