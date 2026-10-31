library(yaml)
config <- read_yaml("config/application_milk_20261031.yml")
source("research/score_validity_20261026/load_joint.R")
output <- "results/preflight/milk_scs_20261031/execution_environment.txt"
lines <- c(paste("Environment recorded", format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE)),
  paste("Execution version", config$execution_version), trimws(capture.output(sessionInfo()), which = "right"))
writeLines(lines, output)
if (!file.copy(output, file.path(config$backup_directory, "summary"), overwrite = TRUE))
  stop("Environment record backup failed.")
