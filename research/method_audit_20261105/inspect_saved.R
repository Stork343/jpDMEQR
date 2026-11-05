library(data.table)
library(jsonlite)
library(yaml)
library(matrixStats)
source("research/clinical_search_20261105/panel.R")
source("research/empirical_selection_20261101/prediction_core.R")
source("research/empirical_selection_20261101/ridge_direction.R")
configurations <- c("config/application_ra_westchina_20261105.yml",
  "config/application_pregnancy_detroit_20261105.yml", "config/application_gnhs_baseline_bmi_20261105.yml")
output <- "results/preflight/method_audit_20261105"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
records <- list()
for (path in configurations) {
  config <- read_yaml(path)
  panel <- cp_panel(config$input_directory)
  for (fold in seq_len(config$outer_folds)) for (method in c("JOINT_DIRECT", "JOINT_OFFSET")) {
    rows <- fread(file.path(config$output_root, paste0("outer_", fold), "training.csv"))
    for (tau in unlist(config$quantiles)) {
      directory <- file.path(config$output_root, paste0("outer_", fold), "models", sprintf("%s_tau%03d", method, as.integer(tau * 100)))
      fit <- readRDS(file.path(directory, "fit.rds"))
      features <- if (method == "JOINT_OFFSET") cp_offset_features(rows, panel, fit$recipe) else cp_features(rows, panel, fit$recipe)
      clinical <- which(startsWith(names(fit$beta), "C_") | startsWith(names(fit$beta), "M_C_"))
      chosen <- which(abs(fit$beta) > 1e-9)
      setting <- as.data.table(fit$grid)[fit$choice]
      support <- unique(c(which(features$penalty_factor == 0),
        vapply(fit$final_diagnostics[[1L]], `[[`, integer(1), "selected_coordinate")))
      stopifnot(all(chosen %in% support), length(support) == setting$size + 1L)
      records[[length(records) + 1L]] <- data.table(study_id = config$study_id, fold, method, tau,
        decision_quantile = tau == config$decision_quantile, training_patients = features$n,
        variables = ncol(features$X), Lambda = setting$Lambda, selected_additions = setting$size,
        maximum_additions = max(fit$grid$size), clinical_columns = length(clinical),
        clinical_forced = sum(features$penalty_factor[clinical] == 0),
        clinical_selected = sum(clinical %in% support),
        proteins_selected = sum(startsWith(names(fit$beta)[support], "P_")),
        missing_indicators_selected = sum(startsWith(names(fit$beta)[support], "M_")),
        direction_selected = sum(startsWith(names(fit$beta)[support], "RIDGE_")),
        omitted_clinical_columns = paste(names(fit$beta)[setdiff(clinical, support)], collapse = ";"),
        calibration = unname(fit$calibration[[1L]][1L]))
    }
  }
}
table <- rbindlist(records)
fwrite(table, file.path(output, "saved_selection_audit.csv"))
summary <- table[decision_quantile == TRUE, .(fits = .N, maximum_grid_Lambda_selected = sum(Lambda == 2),
  minimum_selected_additions = min(selected_additions), median_selected_additions = median(selected_additions),
  maximum_selected_additions = max(selected_additions),
  fits_with_omitted_clinical_columns = sum(clinical_selected < clinical_columns),
  median_selected_proteins = median(proteins_selected)), by = .(study_id, method)]
fwrite(summary, file.path(output, "selection_summary.csv"))
print(summary)
