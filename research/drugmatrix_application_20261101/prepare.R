library(data.table)
library(jsonlite)

directory <- "data/derived/drugmatrix_20261101"
output <- "results/preflight/empirical_selection_20261101/drugmatrix_metadata"
if (file.exists(file.path(directory, "panel.rds"))) stop("数据对象已经存在")
metadata <- fread(file.path(directory, "cohort_metadata.csv"))
stopifnot(!anyDuplicated(metadata$ARRAY_ID), !anyDuplicated(metadata$ANIMAL_ID))
path <- "data/raw/drugmatrix_20261101/extracted/liver/liver for genespring -2s to blank_quantile_genespring.txt"
header <- names(fread(path, nrows = 0L))
stopifnot(header[1L] == "probe", all(metadata$ARRAY_ID %in% header))
expression <- fread(path, select = c("probe", metadata$ARRAY_ID), na.strings = c("", "NA"))
stopifnot(!anyDuplicated(expression$probe))
X <- t(as.matrix(expression[, -1L, with = FALSE]))
dimnames(X) <- list(metadata$ARRAY_ID, expression$probe)
stopifnot(all(is.finite(X) | is.na(X)))
study_date <- as.Date(metadata$STUDY_DATE, format = "%m/%d/%Y")
stopifnot(!anyNA(study_date), all(metadata$DOSE_clinical > 0), all(metadata$DURATION_clinical > 0))
clinical <- data.frame(log_dose = log(metadata$DOSE_clinical),
  log_duration = log(metadata$DURATION_clinical), study_year = factor(format(study_date, "%Y")))
C <- model.matrix(~ log_dose + log_duration + study_year, clinical)
rownames(C) <- metadata$ARRAY_ID
rows <- metadata[, .(sample_id = ARRAY_ID, animal_id = ANIMAL_ID,
                      compound_id, study_id = STUDY_ID, y = log(ALT_numeric),
                      alt_u_l = ALT_numeric, role, exploration_role, inner_fold)]
setorder(rows, compound_id, sample_id)
saveRDS(list(X = X, C = C), file.path(directory, "panel.rds"), compress = FALSE)
saveRDS(rows[role == "discovery"], file.path(directory, "discovery.rds"))
saveRDS(rows[role == "validation"], file.path(directory, "validation_reserved.rds"))
fwrite(rows[, .(animals = .N, compounds = uniqueN(compound_id)), by = .(role, exploration_role)],
       file.path(output, "matrix_cohort.csv"))
write_json(list(status = "prepared", arrays = nrow(X), probes = ncol(X),
                missing_predictor_fraction = mean(is.na(X)), expression_scale = "supplied_normalized_values",
                expression_minimum = min(X, na.rm = TRUE), expression_maximum = max(X, na.rm = TRUE),
                adjustment_columns = colnames(C), maximum_animals_per_compound = max(table(rows$compound_id)),
                actual_response = "natural_log_of_ALT_U_per_L", response_imputed = FALSE,
                formal_validation_used = FALSE, manuscript_adoption = FALSE),
           file.path(output, "preparation.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
cat(nrow(X), "animals", ncol(X), "probes", uniqueN(rows$compound_id), "compounds\n")
