library(dplyr)
root <- "results/raw/simulation/submission_20261027_pilot_joint"
request <- readRDS(file.path(root, "request.rds"))
paths <- list.files(file.path(root, "replications"), "^result[.]rds$", recursive = TRUE, full.names = TRUE)
if (length(paths) != nrow(request$jobs)) stop("Incomplete joint pilot.")
metrics <- coordinates <- vector("list", length(paths))
for (i in seq_along(paths)) {
  result <- readRDS(paths[i])
  metrics[[i]] <- result$metrics
  coordinates[[i]] <- result$coordinates
}
metrics <- bind_rows(metrics)
coordinates <- bind_rows(coordinates)
availability <- metrics %>% group_by(experiment_id, method_id) %>% summarise(
  attempted = n(), converged_fraction = mean(converged), maximum_kkt = if (any(is.finite(kkt_residual))) max(kkt_residual, na.rm = TRUE) else NA_real_,
  median_runtime = median(runtime_sec), failures = sum(!status %in% c("ok", "warning")), .groups = "drop")
replicate <- coordinates %>% mutate(available = is.finite(estimated_se),
  covered_all = available & !is.na(covered) & covered) %>%
  group_by(experiment_id, method_id, coordinate_type, replicate) %>%
  summarise(coverage = mean(covered_all), availability = mean(available), .groups = "drop")
coverage <- replicate %>% group_by(experiment_id, method_id, coordinate_type) %>%
  summarise(mcse = sd(coverage) / sqrt(n()), coverage = mean(coverage),
            availability = mean(availability), replications = n(), .groups = "drop")
output <- "results/tables/submission_20261027/pilot"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
write.csv(availability, file.path(output, "availability.csv"), row.names = FALSE)
write.csv(coverage, file.path(output, "coverage.csv"), row.names = FALSE)
write.csv(coordinates, file.path(output, "coordinate_records.csv"), row.names = FALSE)
gate <- availability %>% filter(method_id == "JOINT-PROFILE")
numeric_pass <- all(gate$converged_fraction >= .98 & gate$maximum_kkt <= 1.1e-8)
jsonlite::write_json(list(numerical_pass = numeric_pass, registered_replications = nrow(request$jobs),
  observed_replications = length(paths), proposed = gate, coverage = filter(coverage, method_id == "JOINT-PROFILE"),
  calibration_action = "none", completed = as.character(Sys.time())), file.path(output, "gate.json"), auto_unbox = TRUE, pretty = TRUE)
print(gate)
print(filter(coverage, method_id == "JOINT-PROFILE"))
if (!numeric_pass) stop("Joint pilot numerical acceptance failed.")
