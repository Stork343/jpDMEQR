library(data.table)
library(jsonlite)

output_dir <- "data/derived/drug_response_20261101"
audit_dir <- "results/preflight/empirical_selection_20261101/drug_response"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
output <- file.path(output_dir, "model_roles.csv")
if (file.exists(output)) stop("样本分组已经存在")
metadata <- fread("data/raw/gdsc_20261101/model_list_20261101.csv.gz")
metadata <- metadata[model_type == "Cell Line"]
stopifnot(!anyDuplicated(metadata$model_id), all(nzchar(metadata$patient_id)),
          !anyNA(metadata$patient_id))
parent <- metadata[nzchar(parent_id)]
parent_patient <- metadata$patient_id[match(parent$parent_id, metadata$model_id)]
if (anyNA(parent_patient) || any(parent$patient_id != parent_patient)) {
  stop("衍生细胞系与来源记录需要进一步核查")
}
patients <- sort(unique(metadata$patient_id))
set.seed(20261101L)
assignment <- data.table(patient_id = patients,
                         role = sample(rep(c("discovery", "validation"),
                                           c(floor(0.7 * length(patients)),
                                             length(patients) - floor(0.7 * length(patients))))))
roles <- merge(metadata[, .(model_id, model_name, patient_id, parent_id, tissue,
                           cancer_type, COSMIC_ID, BROAD_ID, CCLE_ID, RRID)],
               assignment, by = "patient_id", all.x = TRUE, sort = FALSE)
setorder(roles, model_id)
stopifnot(!anyNA(roles$role), roles[, max(uniqueN(role)), by = patient_id][, max(V1)] == 1L)
fwrite(roles, output)
counts <- roles[, .(models = .N, patients = uniqueN(patient_id)), by = role]
fwrite(counts, file.path(audit_dir, "role_counts.csv"))
write_json(list(
  status = "roles_fixed_before_response_analysis",
  seed = 20261101L, grouping = "patient_id",
  parent_consistency_checked = TRUE,
  response_used = FALSE,
  model_metadata_source = "https://cog.sanger.ac.uk/cmp/download/model_list_latest.csv.gz",
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
), file.path(audit_dir, "roles.json"), auto_unbox = TRUE, pretty = TRUE)
print(counts)
