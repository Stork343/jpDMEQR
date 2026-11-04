library(data.table)
library(jsonlite)
library(yaml)
source("research/yeast_stress_application_20261104/features.R")
source("research/yeast_proteome_application_20261105/features.R")
source("research/empirical_selection_20261101/prediction_core.R")
source("research/empirical_selection_20261101/ridge_direction.R")
source("research/empirical_selection_20261101/report_prediction.R")

config <- read_yaml("config/development_ridge_offset_series_20261105.yml")
queue <- fread(file.path(config$output_root, "progress.csv"))
stopifnot(nrow(queue) == 35L, all(queue$status == "completed"))
prior <- rbindlist(lapply(unlist(config$assays), function(assay_name)
  data.table(trait = assay_name, method = c("ENET", "JOINT"), status = "completed",
    directory = file.path("results/raw/application/yeast_proteome_20261105/discovery", assay_name,
                          c("ENET_weighted_cv", "JOINT_prediction_tau025")))))
queue <- rbindlist(list(queue, prior), fill = TRUE)
panel <- proteome_panel(config$input_directory)
offset_features <- function(rows, panel, recipe) ep_offset_features(rows, panel, recipe, proteome_features)
feature_functions <- list(JOINT_OFFSET = offset_features, POOLED_OFFSET = offset_features,
                           PRESS_RIDGE = proteome_features, ENET = proteome_features, JOINT = proteome_features)
prediction_functions <- list(JOINT_OFFSET = ep_predict_offset, POOLED_OFFSET = ep_predict_offset,
                              PRESS_RIDGE = ep_predict, ENET = ep_predict, JOINT = ep_predict)
output <- "results/preflight/empirical_selection_20261101/ridge_offset_series"
comparison <- ep_report_prediction(config, queue, panel, proteome_features, output,
  target_method = "JOINT_OFFSET", stage = "method_development", feature_functions = feature_functions,
  prediction_functions = prediction_functions)
background <- fread(file.path(config$output_root, "recipe_summary.csv"))
runtime <- fread(file.path(output, "settings.csv"))
runtime[, background_and_transform_sec := NA_real_]
for (j in seq_len(nrow(runtime))) {
  assay_name <- runtime$assay[j]
  method_name <- runtime$method[j]
  duration <- if (method_name == "PRESS_RIDGE") {
    background[assay == assay_name & panel == "training", sum(background_and_transform_sec)]
  } else if (method_name %in% c("JOINT_OFFSET", "POOLED_OFFSET")) {
    background[assay == assay_name, sum(background_and_transform_sec)]
  } else 0
  runtime[j, background_and_transform_sec := duration]
}
runtime[, total_training_tuning_sec := runtime_sec + background_and_transform_sec]
fwrite(runtime, file.path(output, "runtime_components.csv"))
screen <- comparison[tau == config$decision_quantile & comparator %in% c("ENET", "PRESS_RIDGE", "POOLED_OFFSET"),
  .(comparators = .N, smallest_relative_gain = min(relative_gain),
     smallest_interval_lower = min(lower), further_comparisons_required = all(relative_gain > 0)), by = assay]
stopifnot(all(screen$comparators == 3L))
fwrite(screen, file.path(output, "screening_decision.csv"))
print(screen)
