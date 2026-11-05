library(data.table)
library(jsonlite)
library(yaml)
library(matrixStats)
source("research/clinical_search_20261105/panel.R")
source("research/empirical_selection_20261101/prediction_core.R")

configurations <- c("config/application_ra_westchina_20261105.yml",
  "config/application_pregnancy_detroit_20261105.yml", "config/application_gnhs_baseline_bmi_20261105.yml")
output <- "results/preflight/method_audit_20261105"
checks <- warning_records <- comparator_records <- list()
for (path in configurations) {
  config <- read_yaml(path)
  panel <- cp_panel(config$input_directory)
  for (fold in seq_len(config$outer_folds)) {
    outer <- file.path(config$output_root, paste0("outer_", fold))
    rows <- fread(file.path(outer, "training.csv"))
    original <- readRDS(file.path(outer, "models",
      sprintf("JOINT_DIRECT_tau%03d", as.integer(config$decision_quantile * 100)), "fit.rds"))
    features <- cp_features(rows, panel, original$recipe)
    for (variant in c("profile_selected", "pinball_selected", "profile_full_clinical", "pinball_full_clinical")) {
      saved <- readRDS(file.path("results/raw/method_audit_20261105", config$study_id,
        sprintf("fold_%d_%s.rds", fold, variant)))
      prediction <- drop(features$X %*% saved$fit$beta) * original$recipe$y_scale + original$recipe$y_center
      messages <- character()
      calibration <- withCallingHandlers(
        ep_calibrate(rows, prediction, config$decision_quantile, "strain_id", "condition"),
        warning = function(value) messages <<- c(messages, conditionMessage(value)))
      stopifnot(max(abs(calibration - saved$calibration)) < 1e-8)
      residual <- rows$y - prediction
      cut <- unname(calibration[1L])
      tolerance <- 1e-8 * max(1, abs(cut))
      below <- sum(features$weight * (residual < cut - tolerance)) / features$n
      at_or_below <- sum(features$weight * (residual <= cut + tolerance)) / features$n
      stopifnot(below <= config$decision_quantile + 1e-10,
        at_or_below >= config$decision_quantile - 1e-10)
      checks[[length(checks) + 1L]] <- data.table(study_id = config$study_id, fold, variant,
        calibration_error = max(abs(calibration - saved$calibration)), below, at_or_below,
        tau = config$decision_quantile, solver_status = saved$fit$solution$solver_status,
        solver_iterations = saved$fit$solution$iterations)
      if (length(messages)) warning_records[[length(warning_records) + 1L]] <- data.table(
        study_id = config$study_id, fold, variant, message = messages)
    }
    enet <- readRDS(file.path(outer, "models", "ENET_all", "fit.rds"))
    selected <- abs(enet$beta) > 1e-9
    setting <- as.data.table(enet$grid)[enet$choice]
    comparator_records[[length(comparator_records) + 1L]] <- data.table(study_id = config$study_id,
      fold, method = "ENET", alpha = setting$alpha, penalty_index = setting$penalty_index,
      selected_columns = sum(selected), selected_proteins = sum(selected & startsWith(names(enet$beta), "P_")))
  }
}
fwrite(rbindlist(checks), file.path(output, "calibration_checks.csv"))
fwrite(rbindlist(warning_records), file.path(output, "calibration_warnings.csv"))
fwrite(rbindlist(comparator_records), file.path(output, "enet_selection.csv"))
writeLines(trimws(capture.output(sessionInfo()), which = "right"), file.path(output, "sessionInfo.txt"))
print(rbindlist(comparator_records))
cat("全部保存模型的校准最优性检查通过\n")
