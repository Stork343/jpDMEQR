library(data.table)
library(jsonlite)
library(GEOquery)

directory <- "data/derived/drugmatrix_20261101"
output <- "results/preflight/empirical_selection_20261101/drugmatrix_metadata"
panel <- readRDS(file.path(directory, "panel.rds"))
source_path <- "data/raw/drugmatrix_20261101/extracted/liver/liver for genespring -2s to blank_quantile_genespring.txt"
sample_ids <- rownames(panel$X)[unique(as.integer(quantile(seq_len(nrow(panel$X)), c(0, 0.25, 0.5, 0.75, 1))))]
reference <- fread(source_path, select = c("probe", sample_ids), na.strings = c("", "NA"))
stopifnot(identical(names(reference)[-1L], sample_ids), identical(reference$probe, colnames(panel$X)))
expected <- t(as.matrix(reference[, -1L, with = FALSE]))
dimnames(expected) <- list(sample_ids, reference$probe)
actual <- panel$X[sample_ids, , drop = FALSE]
stopifnot(identical(is.na(actual), is.na(expected)))
error <- max(abs(actual - expected), na.rm = TRUE)
stopifnot(error < 1e-12)
platform <- getGEO(filename = "data/raw/drugmatrix_20261101/GPL5425_family.soft.gz")
annotation <- as.data.table(Table(platform))
stopifnot(all(c("ID", "TYPE", "GB_ACC", "GENE_SYMBOL", "GENE_TITLE") %in% names(annotation)),
          !anyDuplicated(annotation$ID))
annotation <- annotation[ID %in% colnames(panel$X)]
fwrite(annotation, file.path(directory, "probe_annotation.csv"))
allowed <- annotation[TYPE == "DISCOVERY", ID]
stopifnot(length(allowed) >= 1000L)
write_json(list(status = "passed", source_array_identity_error = error,
                arrays_independently_checked = length(sample_ids), probes_checked_per_array = ncol(panel$X),
                probes_in_platform_annotation = nrow(annotation), biological_probes = length(allowed),
                annotation_source = "GPL5425", platform_title = Meta(platform)$title,
                gene_names_require_current_mapping_for_interpretation = TRUE,
                formal_validation_used = FALSE),
           file.path(output, "input_checks.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
print(annotation[, .N, by = TYPE])
