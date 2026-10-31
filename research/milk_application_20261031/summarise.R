library(data.table)
library(yaml)
config <- read_yaml("config/application_milk_20261031.yml")
registry <- fread(file.path(config$results_directory, "registry.csv"))
run_directory <- file.path(config$results_directory, config$execution_version)
tables <- diagnostics <- selections <- list()
for (i in seq_len(nrow(registry))) {
  job <- registry[i]
  directory <- file.path(run_directory, job$case_id)
  if (!file.exists(file.path(directory, "completion.json"))) stop("The registered empirical case is incomplete: ", job$case_id)
  fit <- readRDS(file.path(directory, "fit.rds"))
  table <- fread(file.path(directory, "inference.csv"))
  if (nrow(table) != 1L || table$case_id != job$case_id || table$n != job$n || table$p != job$p)
    stop("The saved inference does not match its registered case.")
  influence <- fit$details[[1L]]$influence
  if (abs(sqrt(sum(influence^2)) / job$n - table$se) > 1e-12)
    stop("The saved uncertainty differs from the cluster influence calculation.")
  tables[[i]] <- table
  diagnostics[[i]] <- fread(file.path(directory, "diagnostics.csv"), integer64 = "double")
  selections[[i]] <- rbindlist(lapply(1:2, function(fold) {
    current <- fit$prepared$folds[[fold]]
    if (length(intersect(current$training_subjects, current$validation_subjects)))
      stop("A saved selection sample overlaps its estimation sample.")
    graph <- current$graphs[[1L]]
    data.table(case_id = job$case_id, fold = fold, training_cows = length(current$training_subjects),
      estimation_cows = length(current$validation_subjects), outcome_additions = current$initial$size,
      selected_dimension = length(graph$support), curvature_dimension = length(graph$curvature_support),
      design_dimension = length(graph$design_support), design_kkt = graph$design_kkt,
      selected_variables = paste(colnames(current$models[[1L]]$data$X), collapse = ";"))
  }))
}
inference <- rbindlist(tables)
numerical <- rbindlist(diagnostics)
selection <- rbindlist(selections)
inference <- merge(inference, selection[, .(dimension_min = min(selected_dimension),
  dimension_max = max(selected_dimension)), by = case_id], by = "case_id", sort = FALSE)
inference <- inference[match(registry$case_id, case_id)]
output <- "results/preflight/milk_scs_20261031"
fwrite(inference, file.path(output, "inference.csv"))
fwrite(numerical, file.path(output, "fit_diagnostics.csv"))
fwrite(selection, file.path(output, "selected_variables.csv"))
manifest <- list(cases_registered = nrow(registry), cases_completed = nrow(inference),
  coefficient_systems = nrow(numerical), maximum_kkt = max(numerical$kkt),
  maximum_nuisance_gradient = max(numerical$nuisance_gradient),
  maximum_prediction_identity_error = max(numerical$prediction_identity_error),
  maximum_hessian_condition = max(numerical$hessian_condition),
  maximum_basis_condition = max(numerical$basis_condition),
  minimum_selected_dimension = min(selection$selected_dimension),
  maximum_selected_dimension = max(selection$selected_dimension),
  median_runtime_seconds = median(inference$runtime_seconds),
  maximum_runtime_seconds = max(inference$runtime_seconds),
  execution_version = config$execution_version,
  summarised = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE))
jsonlite::write_json(manifest, file.path(output, "completion_manifest.json"),
  auto_unbox = TRUE, pretty = TRUE, digits = NA)
backup <- file.path(config$backup_directory, "summary")
dir.create(backup, recursive = TRUE, showWarnings = FALSE)
if (!all(file.copy(list.files(output, full.names = TRUE), backup, overwrite = TRUE))) stop("Summary backup failed.")
print(inference)
print(manifest)
