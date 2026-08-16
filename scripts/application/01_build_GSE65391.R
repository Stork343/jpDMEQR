#!/usr/bin/env Rscript
full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "scripts/application/01_build_GSE65391.R"
root <- if (file.exists("R/application_v2.R")) "." else normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", envir = environment())
source_v2_module(root, "application_v2", envir = environment())

if (!requireNamespace("GEOquery", quietly = TRUE) || !requireNamespace("Biobase", quietly = TRUE)) {
  stop("Install GEOquery and Biobase first.")
}
raw_dir <- file.path(root, "data", "raw", "GSE65391")
interim_dir <- file.path(root, "data", "interim", "GSE65391")
dir.create(interim_dir, recursive = TRUE, showWarnings = FALSE)
matrix_file <- file.path(raw_dir, "GSE65391_series_matrix.txt.gz")
if (!file.exists(matrix_file)) stop("Series Matrix missing. Run 00_download_GSE65391.R first.")

message("Reading Series Matrix; this can take several minutes.")
geo <- GEOquery::getGEO(filename = matrix_file, GSEMatrix = TRUE, getGPL = FALSE)
eset <- if (inherits(geo, "ExpressionSet")) geo else geo[[1]]
obj <- build_gse65391_objects_v2(eset)

saveRDS(eset, file.path(interim_dir, "gse65391_eset.rds"), compress = "xz")
saveRDS(obj$expression, file.path(interim_dir, "expression_processed_probe.rds"), compress = "xz")
write.csv(obj$metadata, file.path(interim_dir, "phenotype_parsed.csv"), row.names = FALSE)
write.csv(obj$duplicate_characteristic_keys %||% data.frame(),
          file.path(interim_dir, "duplicate_characteristic_keys.csv"), row.names = FALSE)

# Data dictionary before any model-specific recoding.
dictionary <- data.frame(
  variable = names(obj$metadata),
  class = vapply(obj$metadata, function(x) paste(class(x), collapse = "/"), character(1)),
  missing = vapply(obj$metadata, function(x) sum(is.na(x)), integer(1)),
  distinct = vapply(obj$metadata, function(x) length(unique(x[!is.na(x)])), integer(1)),
  stringsAsFactors = FALSE
)
write.csv(dictionary, file.path(interim_dir, "metadata_dictionary.csv"), row.names = FALSE)

# Probe annotation. The exact annotation package version is recorded in sessionInfo.
probe_ids <- rownames(obj$expression)
annotation <- data.frame(probe_id = probe_ids, gene_symbol = NA_character_,
                         gene_name = NA_character_, entrez_id = NA_character_,
                         is_control = FALSE, is_multimapping = FALSE,
                         stringsAsFactors = FALSE)
if (requireNamespace("illuminaHumanv4.db", quietly = TRUE) &&
    requireNamespace("AnnotationDbi", quietly = TRUE)) {
  ann_raw <- AnnotationDbi::select(
    illuminaHumanv4.db::illuminaHumanv4.db,
    keys = probe_ids,
    keytype = "PROBEID",
    columns = c("SYMBOL", "GENENAME", "ENTREZID")
  )
  symbol_count <- tapply(ann_raw$SYMBOL, ann_raw$PROBEID,
                         function(x) length(unique(x[!is.na(x) & nzchar(x)])))
  first_nonmissing <- function(x) {
    x <- x[!is.na(x) & nzchar(x)]
    if (length(x)) x[1] else NA_character_
  }
  by_probe <- split(ann_raw, ann_raw$PROBEID)
  ann_one <- lapply(by_probe, function(d) data.frame(
    probe_id = d$PROBEID[1],
    gene_symbol = first_nonmissing(d$SYMBOL),
    gene_name = first_nonmissing(d$GENENAME),
    entrez_id = first_nonmissing(as.character(d$ENTREZID)),
    is_control = FALSE,
    is_multimapping = (symbol_count[d$PROBEID[1]] %||% 0) > 1,
    stringsAsFactors = FALSE
  ))
  ann_one <- do.call(rbind, ann_one)
  annotation <- merge(data.frame(probe_id = probe_ids, stringsAsFactors = FALSE),
                      ann_one, by = "probe_id", all.x = TRUE, sort = FALSE)
  annotation$is_control[is.na(annotation$is_control)] <- FALSE
  annotation$is_multimapping[is.na(annotation$is_multimapping)] <- FALSE
}
write.csv(annotation, file.path(interim_dir, "probe_annotation.csv"), row.names = FALSE)

capture.output(sessionInfo(), file = file.path(interim_dir, "build_sessionInfo.txt"))
message("Built interim objects under: ", interim_dir)
