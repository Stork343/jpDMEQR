library(yaml)
library(jsonlite)
config <- read_yaml("config/application_soynam_20261031.yml")
for (case_id in c("pooled_quantile_lasso_q50", "pooled_quantile_lasso_q25")) {
  root <- file.path(config$results_directory, "prediction", config$execution_version)
  directory <- file.path(root, case_id)
  if (!dir.exists(directory) || file.exists(file.path(directory, "completion.json")) ||
      file.exists(file.path(directory, "predictions.csv"))) stop("The pooled attempt cannot be archived by this numerical revision.")
  write_json(list(status = "superseded_optimizer", test_predictions_calculated = FALSE,
    completed_training_folds = 0L, replacement = "certified_unsmoothed_quantile_lasso_lp",
    archived = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE)),
    file.path(directory, "interruption.json"), auto_unbox = TRUE, pretty = TRUE)
  destination <- file.path(root, "numerical_attempts", paste0(case_id, "_native_hqreg"))
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  if (dir.exists(destination) || !file.rename(directory, destination)) stop("The pooled attempt archive failed.")
  backup <- file.path(config$backup_directory, "prediction", config$execution_version, case_id)
  backup_destination <- file.path(config$backup_directory, "prediction", config$execution_version,
    "numerical_attempts", paste0(case_id, "_native_hqreg"))
  dir.create(backup_destination, recursive = TRUE, showWarnings = FALSE)
  if (!all(file.copy(list.files(destination, full.names = TRUE), backup_destination, overwrite = FALSE)))
    stop("The pooled attempt backup failed.")
  if (dir.exists(backup) && length(list.files(backup)) == 0L) unlink(backup, recursive = TRUE)
}
