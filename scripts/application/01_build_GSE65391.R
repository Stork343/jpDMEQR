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
protected <- file.path(interim_dir, c("gse65391_eset.rds", "expression_processed_probe.rds", "phenotype_parsed.csv"))
if (any(file.exists(protected))) stop("INTERIM_OUTPUT_EXISTS: preserve the existing GSE65391 objects.")

message("Reading Series Matrix; this can take several minutes.")
geo <- GEOquery::getGEO(filename = matrix_file, GSEMatrix = TRUE, getGPL = FALSE)
eset <- if (inherits(geo, "ExpressionSet")) geo else geo[[1]]
obj <- build_gse65391_objects_v2(eset)

# 使用有版本记录的官方平台注释。
annotation_path <- file.path(interim_dir, "probe_annotation_20261022.csv")
if (file.exists(annotation_path)) {
  annotation <- read.csv(annotation_path, stringsAsFactors = FALSE)
  if (anyDuplicated(annotation$probe_id) || !all(rownames(obj$expression) %in% annotation$probe_id) ||
      !all(c("gene_symbol", "mapping_status", "annotation_source") %in% names(annotation))) {
    stop("Official versioned annotation does not cover the observed expression probes.")
  }
} else {
  platform_dir <- file.path(root, "data", "raw", "public_modules", "banchereau2016_20261022")
  curated_file <- file.path(platform_dir, "GPL10558.annot.gz")
  submitter_file <- file.path(platform_dir, "GPL10558_submitter_complete.soft")
  if (!file.exists(curated_file) || !file.exists(submitter_file) ||
      system2("gzip", c("-t", shQuote(curated_file))) != 0L) {
    stop("Complete official GPL10558 annotation sources are required.")
  }
  curated <- GEOquery::getGEO(filename = curated_file, AnnotGPL = TRUE, destdir = platform_dir)
  vendor <- GEOquery::getGEO(filename = submitter_file, destdir = platform_dir)
  annotation <- assemble_gpl10558_probe_annotation_v2(
    rownames(obj$expression), GEOquery::Table(curated), GEOquery::Table(vendor),
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10558/annot/GPL10558.annot.gz",
    as.character(GEOquery::Meta(curated)$date)
  )
  write.csv(annotation, annotation_path, row.names = FALSE)
}

saveRDS(eset, file.path(interim_dir, "gse65391_eset.rds"), compress = "xz")
saveRDS(obj$expression, file.path(interim_dir, "expression_processed_probe.rds"), compress = "xz")
write.csv(obj$metadata, file.path(interim_dir, "phenotype_parsed.csv"), row.names = FALSE)
write.csv(obj$duplicate_characteristic_keys %||% data.frame(),
          file.path(interim_dir, "duplicate_characteristic_keys.csv"), row.names = FALSE)

# 保存metadata字段字典。
dictionary <- data.frame(
  variable = names(obj$metadata),
  class = vapply(obj$metadata, function(x) paste(class(x), collapse = "/"), character(1)),
  missing = vapply(obj$metadata, function(x) sum(is.na(x)), integer(1)),
  distinct = vapply(obj$metadata, function(x) length(unique(x[!is.na(x)])), integer(1)),
  stringsAsFactors = FALSE
)
write.csv(dictionary, file.path(interim_dir, "metadata_dictionary.csv"), row.names = FALSE)

capture.output(sessionInfo(), file = file.path(interim_dir, "build_sessionInfo.txt"))
message("Built interim objects under: ", interim_dir)
