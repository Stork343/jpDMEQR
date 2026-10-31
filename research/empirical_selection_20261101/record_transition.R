library(data.table)
library(jsonlite)
root <- "results/raw/application/soynam_yield_20261031"
prediction <- file.path(root, "prediction/paired_qp_20261031")
status <- fread(file.path(prediction, "status.csv"))
status[state == "running", `:=`(state = "terminated_not_adopted", exit_status = 143L)]
destination <- "results/preflight/empirical_selection_20261101"
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
fwrite(status, file.path(destination, "previous_application_execution.csv"))
record <- list(application_id = "soynam_yield_20261031", manuscript_adoption = FALSE,
  reason = "User requires demonstrated comparative advantage or a verified substantive new finding.",
  existing_results_retained = TRUE, automatic_manuscript_updates = "terminated",
  recorded = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE))
write_json(record, file.path(destination, "selection_transition.json"), auto_unbox = TRUE, pretty = TRUE)
write_json(record, file.path(root, "selection_status.json"), auto_unbox = TRUE, pretty = TRUE)
print(status[, .N, by = state])
