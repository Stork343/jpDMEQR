library(yaml)
library(jsonlite)
config <- read_yaml("config/application_soynam_20261031.yml")
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S", tz = "Asia/Shanghai")
destination <- file.path(config$backup_directory, "deliveries", stamp)
if (dir.exists(destination)) stop("The delivery snapshot already exists.")
dir.create(destination, recursive = TRUE)
files <- c("manuscript/sjs_profile_quantile_main.pdf", "manuscript/sjs_profile_quantile_supplement.pdf",
  "research/soynam_application_20261031/RESULTS_zh.md", "config/application_soynam_20261031.yml",
  list.files("manuscript/generated/soynam_profile", full.names = TRUE),
  list.files("results/tables/soynam_20261031", full.names = TRUE),
  list.files("results/preflight/soynam_20261031/joint_prediction", full.names = TRUE))
if (!all(file.exists(files))) stop("A requested delivery artifact is missing.")
for (file in files) {
  target <- file.path(destination, file)
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(file, target, overwrite = FALSE) || file.info(file)$size != file.info(target)$size)
    stop("The delivery artifact copy failed.")
}
write_json(list(created = as.character(Sys.time()), files = files,
  source_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE)),
  file.path(destination, "snapshot.json"), auto_unbox = TRUE, pretty = TRUE)
cat(destination, "\n")
