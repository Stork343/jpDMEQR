library(data.table)
library(jsonlite)

directory <- "results/raw/application/ridge_offset_series_20261105/NaCl/JOINT_OFFSET_tau025"
status <- read_json(file.path(directory, "status.json"), simplifyVector = TRUE)
stopifnot(status$status == "completed")
previous <- fread("results/raw/application/calibrated_profile_20261105/OFFSET_CALIBRATED/predictions.csv")
current <- fread(file.path(directory, "predictions.csv"))
stopifnot(identical(previous$sample_id, current$sample_id))
error <- max(abs(previous$prediction - current$prediction), abs(previous$loss - current$loss))
stopifnot(is.finite(error), error < 1e-9)
fit <- readRDS(file.path(directory, "fit.rds"))
stopifnot(fit$calibration_mode == "training")
output <- "results/preflight/empirical_selection_20261101/ridge_direction"
write_json(list(status = "passed", compared_actual_predictions = nrow(current), maximum_error = error,
                formal_validation_used = FALSE), file.path(output, "calibrated_selection_equivalence.json"),
           auto_unbox = TRUE, pretty = TRUE, digits = NA)
cat("完整预测流程与保存路径计算结果一致\n")
