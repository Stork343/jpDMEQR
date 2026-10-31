library(data.table)
library(jsonlite)

source_file <- "data/raw/idrr_20261101/iDRR_data.rds"
output_dir <- "results/preflight/empirical_selection_20261101/idrr"
interim_dir <- "data/interim/drug_response_20261101"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(interim_dir, recursive = TRUE, showWarnings = FALSE)
if (file.exists(file.path(output_dir, "structure.json"))) stop("数据核查记录已经存在")
source <- readRDS(source_file)
response <- as.data.table(source$dose_response)
drugs <- as.data.table(source$drugs)
fwrite(drugs, file.path(output_dir, "compound_metadata.csv"))
roles <- fread("data/derived/drug_response_20261101/model_roles.csv")
mapping <- unique(roles[nzchar(RRID), .(cello_id = RRID, patient_id, role)])
ambiguities <- mapping[, .N, by = cello_id][N > 1L]
fwrite(ambiguities, file.path(output_dir, "identifier_ambiguities.csv"))
if (nrow(ambiguities)) stop("Cellosaurus 标识存在多个患者来源")
availability <- response[, .(rows = .N, cells = uniqueN(cello_id), compounds = uniqueN(cid)), by = dataset]
setorder(availability, -rows)
fwrite(availability, file.path(output_dir, "study_availability.csv"))
selected_drugs <- drugs[grepl("dabrafenib|plx.?4720|az.?628|sb.?590885|hg.?6.?64.?1", drug,
                              ignore.case = TRUE)]
if (!nrow(selected_drugs)) stop("药物名称需要进一步核查")
fwrite(selected_drugs, file.path(output_dir, "selected_compound_metadata.csv"))
selected <- response[cid %in% selected_drugs$cid]
selected <- merge(selected, selected_drugs[, .(cid, drug)], by = "cid", all.x = TRUE)
selected <- merge(selected, mapping, by = "cello_id", all.x = TRUE, sort = FALSE)
counts <- selected[, .(rows = .N, cells = uniqueN(cello_id),
                        patients = uniqueN(patient_id[!is.na(patient_id)]),
                        matched_rows = sum(!is.na(patient_id)),
                        discovery_rows = sum(role == "discovery", na.rm = TRUE),
                        validation_rows = sum(role == "validation", na.rm = TRUE)),
                    by = .(dataset, drug)]
fwrite(counts, file.path(output_dir, "candidate_availability.csv"))
discovery <- selected[role == "discovery"]
quality <- discovery[, .(
  rows = .N, patients = uniqueN(patient_id), cells = uniqueN(cello_id),
  missing_response = sum(!is.finite(viability)),
  minimum_response = min(viability, na.rm = TRUE),
  maximum_response = max(viability, na.rm = TRUE),
  response_zero = mean(viability == 0, na.rm = TRUE),
  response_hundred = mean(viability == 100, na.rm = TRUE),
  response_negative = mean(viability < 0, na.rm = TRUE),
  response_above_hundred = mean(viability > 100, na.rm = TRUE),
  missing_dose = sum(!is.finite(dose)), nonpositive_dose = sum(dose <= 0, na.rm = TRUE)
), by = .(dataset, drug)]
fwrite(quality, file.path(output_dir, "discovery_measurement_quality.csv"))
repetitions <- discovery[, .(rows = .N, doses = uniqueN(dose)), by = .(dataset, drug, patient_id)]
structure <- repetitions[, .(patients = .N, minimum_rows = min(rows), median_rows = median(rows),
                               maximum_rows = max(rows), minimum_doses = min(doses),
                               median_doses = median(doses), maximum_doses = max(doses)),
                          by = .(dataset, drug)]
fwrite(structure, file.path(output_dir, "discovery_cluster_structure.csv"))
saveRDS(discovery, file.path(interim_dir, "idrr_braf_discovery.rds"))
saveRDS(selected[role == "validation"], file.path(interim_dir, "idrr_braf_validation_reserved.rds"))
write_json(list(status = "structure_audited", source_file = source_file,
                source_rows = nrow(response), candidate_rows = nrow(selected),
                response_association_fitted = FALSE, validation_results_used = FALSE,
                manuscript_adoption = FALSE, description = source$description),
           file.path(output_dir, "structure.json"), auto_unbox = TRUE, pretty = TRUE)
print(quality)
