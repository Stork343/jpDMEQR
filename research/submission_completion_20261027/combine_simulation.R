library(data.table)
base <- "results/tables/submission_20261028/collections"
families <- c("submission_20261027_main_joint", paste0("submission_20261028_main_",
  c("pooled", "mixed", "qgee", "bias_adjusted", "sqr")))
paths <- file.path(base, families, "records.rds")
if (any(!file.exists(paths))) stop("The complete simulation collections are required.")
collections <- lapply(paths, readRDS)
metrics <- rbindlist(lapply(collections, `[[`, "metrics"), fill = TRUE)
coordinates <- rbindlist(lapply(collections, `[[`, "coordinates"), fill = TRUE)
theory <- rbindlist(lapply(collections, `[[`, "theory"), fill = TRUE)
specification <- yaml::read_yaml("config/joint_profile_submission_20261027.yml")
registry <- fread(specification$generating_registry)[module != "G_computation"]
mapping <- unlist(specification$method_mapping)
registry[, methods := vapply(strsplit(methods, "|", fixed = TRUE), function(values) {
  matched <- match(values, names(mapping))
  values[!is.na(matched)] <- mapping[matched[!is.na(matched)]]
  paste(values, collapse = "|")
}, character(1))]
for (id in names(specification$target_mode_overrides)) registry[experiment_id == id,
  target_mode := specification$target_mode_overrides[[id]]]
registry[, budget_multiplier := 1]
for (id in names(specification$budget_multipliers)) registry[experiment_id == id,
  budget_multiplier := specification$budget_multipliers[[id]]]
expected <- rbindlist(lapply(seq_len(nrow(registry)), function(i) {
  row <- registry[i]
  CJ(experiment_id = row$experiment_id, replicate = seq_len(row$replications),
    method_id = strsplit(row$methods, "|", fixed = TRUE)[[1L]])
}))
keys <- c("experiment_id", "replicate", "method_id")
if (anyDuplicated(metrics[, ..keys]) || nrow(metrics) != nrow(expected) ||
    nrow(fsetdiff(expected, metrics[, ..keys])) || nrow(fsetdiff(metrics[, ..keys], expected)))
  stop("The final statistical study differs from its complete registry.")
if (any(metrics[, .(seeds = uniqueN(seed)), by = .(experiment_id, replicate)]$seeds != 1L))
  stop("A paired simulation sample has different generating seeds.")
inference <- unique(expected$method_id[startsWith(expected$method_id, "JOINT-PROFILE") |
  expected$method_id %in% c("SQR-DEBIASED-IID", "BIAS-ADJ-LQMM")])
counts <- coordinates[method_id %in% inference, .(coordinates = .N), by = keys]
counts <- merge(counts, registry[, .(experiment_id, target_coordinate_count)], by = "experiment_id")
if (nrow(counts) != nrow(expected[method_id %in% inference]) || any(counts$coordinates != counts$target_coordinate_count))
  stop("Requested simulation coordinates are incomplete.")
output <- file.path(base, "statistical_study")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
saveRDS(list(metrics = metrics, coordinates = coordinates, theory = theory, registered = expected,
  registry = registry, source = paths,
  implementation_commit = setNames(lapply(collections, `[[`, "implementation_commit"), families)),
  file.path(output, "records.rds"))
fwrite(registry, file.path(output, "config_snapshot.csv"))
accounting <- rbindlist(lapply(file.path(base, families, "accounting.csv"), fread))
fwrite(accounting, file.path(output, "accounting.csv"))
status_counts <- accounting[, lapply(.SD, sum), by = method_id,
  .SDcols = c("registered", "converged_count", "numerical_failures", "component_nonconvergence", "domain_unavailable", "unavailable")]
fwrite(status_counts, file.path(output, "method_accounting.csv"))
jsonlite::write_json(list(settings = nrow(registry), datasets = sum(registry$replications),
  registered_method_records = nrow(expected), coordinate_records = nrow(coordinates),
  source_collections = paths, generated = as.character(Sys.time())),
  file.path(output, "manifest.json"), auto_unbox = TRUE, pretty = TRUE)
print(status_counts)
