args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(".", mustWork = TRUE)
run_id <- if (length(args)) args[[1L]] else stop("An immutable check run ID is required.")
out <- file.path(root, "results", "scratch", "source_registry_checks_20261022", run_id)
if (dir.exists(out)) stop("The check directory already exists.")
dir.create(out, recursive = TRUE)
files <- unlist(lapply(c("R", "scripts", "tests"), function(directory) {
  list.files(file.path(root, directory), pattern = "\\.[Rr]$", recursive = TRUE, full.names = TRUE)
}))
parsed <- lapply(files, function(path) {
  result <- parse(file = path, keep.source = FALSE)
  data.frame(file = substring(path, nchar(root) + 2L), expressions = length(result), pass = TRUE)
})
write.csv(do.call(rbind, parsed), file.path(out, "parsed_sources.csv"), row.names = FALSE)
source(file.path(root, "scripts", "00_source_v2.R"))
for (module in c("profile_v2", "simulation_v2", "metrics_v2", "benchmark_adapters_v2", "application_v2")) {
  source_v2_module(root, module, environment())
}
source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"))
checks <- list(
  list("simulation_pair_main_20261022.yml", "main", 78L),
  list("simulation_pair_pilot_20261022.yml", "pilot", 6L),
  list("simulation_pair_preflight_20261022.yml", "preflight", 10L),
  list("simulation_application_method_preflight_pair_20261022.csv", "preflight", 2L)
)
records <- lapply(checks, function(entry) {
  path <- file.path(root, "config", entry[[1L]])
  registry <- read_simulation_registry_v2(path)
  check <- validate_registry_contract_v2(root, path, entry[[2L]])
  if (!check$pass) stop(paste(check$problems, collapse = "; "))
  if (nrow(registry) != entry[[3L]] || anyDuplicated(registry$experiment_id)) {
    stop("The registry row count or identifiers changed: ", path)
  }
  data.frame(file = entry[[1L]], phase = entry[[2L]], scenarios = nrow(registry),
             requested_replications = sum(registry$replications), pass = check$pass)
})
table <- do.call(rbind, records)
write.csv(table, file.path(out, "registry_checks.csv"), row.names = FALSE)
manifest <- list(pass = all(table$pass), commit_sha = current_commit_v2(root),
                 parsed_files = length(files), registry_checks = table,
                 completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
write_json_v2(manifest, file.path(out, "manifest.json"))
print(manifest)
