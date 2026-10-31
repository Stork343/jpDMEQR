library(data.table)
library(jsonlite)
library(gdscIC50)

root <- "data/raw/gdsc_20261101"
output_dir <- "data/derived/drug_response_20261101"
audit_dir <- "results/preflight/empirical_selection_20261101/gdsc1"
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
if (file.exists(file.path(output_dir, "gdsc1_braf_discovery.rds"))) stop("数据对象已经存在")
roles <- fread(file.path(output_dir, "model_roles.csv"))
compounds <- fread(file.path(root, "screened_compounds_rel_8.5.csv"))
compound_ids <- c(29L, 159L, 1036L, 1061L, 1371L, 1373L)
stopifnot(all(compound_ids %in% compounds$DRUG_ID))

expression_file <- file.path(root, "extracted/rnaseq_fpkm_20191101.csv")
header <- fread(expression_file, header = FALSE, nrows = 4L)
expression <- fread(expression_file, header = FALSE, skip = 4L)
stopifnot(identical(as.character(header[[1L]]), c("model_id", "model_name", "dataset_name", "gene_id")))
model_ids <- as.character(unlist(header[1L, 3:ncol(header), with = FALSE], use.names = FALSE))
expression_source <- as.character(unlist(header[3L, 3:ncol(header), with = FALSE], use.names = FALSE))
stopifnot(!anyDuplicated(model_ids), !anyDuplicated(expression[[1L]]))
annotation <- fread(file.path(root, "gene_identifiers_20261101.csv.gz"))
annotation <- annotation[nzchar(uniprot_id) & nzchar(hgnc_symbol)]
gene_rows <- which(expression[[1L]] %in% annotation$gene_id)
gene_ids <- expression[[1L]][gene_rows]
gene_info <- annotation[match(gene_ids, gene_id)]
stopifnot(identical(gene_info$gene_id, gene_ids), !anyDuplicated(gene_ids))
values <- t(as.matrix(expression[gene_rows, 3:ncol(expression), with = FALSE]))
dimnames(values) <- list(model_ids, gene_ids)
stopifnot(all(values[!is.na(values)] >= 0))
expression_metadata <- data.table(model_id = model_ids, expression_source = expression_source)
saveRDS(list(fpkm = values, genes = gene_info, models = expression_metadata),
        file.path(output_dir, "expression.rds"), compress = FALSE)
rm(expression, values)
gc()

raw <- fread(file.path(root, "extracted/GDSC1_public_raw_data_27Oct23/GDSC1_public_raw_data_27Oct23.csv"))
source_rows <- nrow(raw)
selected_scans <- unique(raw[DRUG_ID %in% compound_ids, SCAN_ID])
raw <- raw[SCAN_ID %in% selected_scans & (DRUG_ID %in% compound_ids | TAG %in% c("NC-0", "B", "FAIL"))]
plate_metadata <- unique(raw[, .(SCAN_ID, SANGER_MODEL_ID, ASSAY, DURATION, SEEDING_DENSITY)])
stopifnot(!anyDuplicated(plate_metadata$SCAN_ID))
initial_rows <- nrow(raw)
raw <- as.data.table(removeMissingDrugs(removeFailedDrugs(as.data.frame(raw))))
quality_removed_rows <- initial_rows - nrow(raw)
controls <- raw[TAG %in% c("NC-0", "B"), .(count = .N, mean = mean(INTENSITY)), by = .(SCAN_ID, TAG)]
fwrite(controls[, .(plates = .N, minimum_count = min(count), median_count = median(count)), by = TAG],
       file.path(audit_dir, "control_structure.csv"))
control_wide <- dcast(controls, SCAN_ID ~ TAG, value.var = "mean")
stopifnot(all(is.finite(control_wide$B)), all(is.finite(control_wide[["NC-0"]])),
          all(control_wide[["NC-0"]] > control_wide$B))
normalized <- as.data.table(normalizeData(as.data.frame(raw), trim = FALSE,
                                        neg_control = "NC-0", pos_control = "B"))
normalized <- normalized[DRUG_ID_lib %in% compound_ids]
stopifnot(all(normalized$treatment == "S"), !anyDuplicated(normalized[, .(SCAN_ID, POSITION)]),
          all(is.finite(normalized$normalized_intensity)), all(is.finite(normalized$CONC)),
          all(normalized$CONC > 0))
normalized <- merge(normalized, plate_metadata, by = "SCAN_ID", all.x = TRUE)
normalized <- merge(normalized, roles, by.x = "SANGER_MODEL_ID", by.y = "model_id", all.x = TRUE)
normalized <- merge(normalized, expression_metadata, by.x = "SANGER_MODEL_ID", by.y = "model_id", all.x = TRUE)
normalized <- merge(normalized, compounds[, .(DRUG_ID_lib = DRUG_ID, DRUG_NAME, SCREENING_SITE)],
                    by = "DRUG_ID_lib", all.x = TRUE)
availability <- normalized[, .(rows = .N, models = uniqueN(SANGER_MODEL_ID)),
                           by = .(DRUG_NAME, role, expression_available = !is.na(expression_source))]
fwrite(availability, file.path(audit_dir, "cohort_availability.csv"))
eligible <- normalized[!is.na(role) & !is.na(expression_source)]
stopifnot(nrow(eligible) > 0L)
setorder(eligible, DRUG_NAME, patient_id, SANGER_MODEL_ID, SCAN_ID, lib_drug, CONC, POSITION)
discovery <- eligible[role == "discovery"]
validation <- eligible[role == "validation"]
stopifnot(!length(intersect(discovery$patient_id, validation$patient_id)))
saveRDS(discovery, file.path(output_dir, "gdsc1_braf_discovery.rds"))
saveRDS(validation, file.path(output_dir, "gdsc1_braf_validation_reserved.rds"))
counts <- eligible[, .(rows = .N, patients = uniqueN(patient_id), models = uniqueN(SANGER_MODEL_ID),
                       experiments = uniqueN(paste(SCAN_ID, lib_drug))), by = .(DRUG_NAME, role)]
fwrite(counts, file.path(audit_dir, "cohort_counts.csv"))
quality <- discovery[, .(
  rows = .N, patients = uniqueN(patient_id),
  minimum = min(normalized_intensity), q25 = quantile(normalized_intensity, 0.25),
  median = median(normalized_intensity), q75 = quantile(normalized_intensity, 0.75),
  maximum = max(normalized_intensity),
  exactly_zero = mean(normalized_intensity == 0), exactly_one = mean(normalized_intensity == 1),
  below_zero = mean(normalized_intensity < 0), above_one = mean(normalized_intensity > 1)
), by = DRUG_NAME]
fwrite(quality, file.path(audit_dir, "discovery_response_quality.csv"))
write_json(list(
  status = "prepared", source_rows = source_rows, quality_removed_rows = quality_removed_rows,
  genes_with_protein_annotation = length(gene_ids), expression_models = length(model_ids),
  preprocessing = "gdscIC50_0.99.4_NC0_B_trim_FALSE",
  validation_response_associations_used = FALSE, manuscript_adoption = FALSE
), file.path(audit_dir, "preparation.json"), pretty = TRUE, auto_unbox = TRUE)
print(counts)
print(quality)
