library(data.table)
library(jsonlite)
library(yaml)

root <- "data/derived/drug_response_20261101"
output <- file.path(root, "discovery_cases.rds")
if (file.exists(output)) stop("探索数据已经准备完成")
config <- read_yaml("config/discovery_gdsc_20261101.yml")
roles <- fread(file.path(root, "model_roles.csv"))
patients <- sort(unique(roles[role == "discovery", patient_id]))
set.seed(config$exploration_seed)
allocation <- data.table(patient_id = patients,
  exploration_role = sample(rep(c("training", "evaluation"),
    c(floor(config$exploration_training_fraction * length(patients)),
      length(patients) - floor(config$exploration_training_fraction * length(patients))))))
allocation[exploration_role == "training", inner_fold := sample(rep(seq_len(config$inner_folds), length.out = .N))]
rows <- readRDS(file.path(root, "gdsc1_braf_discovery.rds"))
stopifnot(all(rows$role == "discovery"))
curves <- rows[, .(doses = uniqueN(CONC), rows = .N,
                   date = min(DATE_CREATED),
                   complete = setequal(dose, paste0("D", 1:5)) || setequal(dose, paste0("D", 1:9))),
               by = .(DRUG_NAME, patient_id, SANGER_MODEL_ID, SCAN_ID, lib_drug)]
setorder(curves, DRUG_NAME, patient_id, date, SCAN_ID, lib_drug)
chosen <- curves[complete == TRUE, .SD[1L], by = .(DRUG_NAME, patient_id)]
chosen[, analysis_curve := TRUE]
rows <- merge(rows, chosen[, .(DRUG_NAME, patient_id, SCAN_ID, lib_drug, analysis_curve)],
              by = c("DRUG_NAME", "patient_id", "SCAN_ID", "lib_drug"), all.x = TRUE, sort = FALSE)
rows[is.na(analysis_curve), analysis_curve := FALSE]
rows <- merge(rows, allocation, by = "patient_id", all.x = TRUE)
stopifnot(!anyNA(rows$exploration_role),
          !length(intersect(rows[exploration_role == "training", patient_id],
                            rows[exploration_role == "evaluation", patient_id])))
counts <- rows[, .(rows = .N, patients = uniqueN(patient_id)),
               by = .(DRUG_NAME, exploration_role, analysis_curve)]
audit_dir <- "results/preflight/empirical_selection_20261101/gdsc_discovery"
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(counts, file.path(audit_dir, "sample_counts.csv"))
fwrite(curves[, .(experiments = .N), by = .(DRUG_NAME, complete, doses)],
       file.path(audit_dir, "curve_completeness.csv"))
fwrite(allocation, file.path(root, "exploration_roles.csv"))
saveRDS(rows, output)
write_json(list(status = "prepared", response_associations_used = FALSE,
                formal_validation_used = FALSE, seed = config$exploration_seed,
                curve_rule = config$curve, remaining_experiments = config$remaining_experiments),
           file.path(audit_dir, "preparation.json"), auto_unbox = TRUE, pretty = TRUE)
print(counts[analysis_curve == TRUE])
