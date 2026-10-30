library(data.table)
library(Biobase)
release <- commandArgs(trailingOnly = TRUE)[1L]
path <- file.path("data/raw/GSE65391", paste0("GSE65391_non-normalized_data_Illumina_HT12_V4_", release, ".txt.gz"))
raw <- fread(path, check.names = FALSE, nThread = 1L, showProgress = FALSE)
eset <- readRDS("data/interim/GSE65391/gse65391_eset.rds")
phenotype <- pData(eset)
arrays <- trimws(sub("^Sample name: ", "", phenotype$description))
signals <- seq.int(2L, ncol(raw), 2L)
print(list(first_column = names(raw)[1L], rows = nrow(raw), columns = ncol(raw),
  duplicate_probes = anyDuplicated(raw[[1L]]), detection_headers = unique(names(raw)[signals + 1L]),
  observed_probes = nrow(exprs(eset)), raw_only = length(setdiff(raw[[1L]], rownames(exprs(eset)))),
  processed_only = length(setdiff(rownames(exprs(eset)), raw[[1L]])),
  unmatched_arrays = setdiff(names(raw)[signals], arrays)))
print(head(setdiff(raw[[1L]], rownames(exprs(eset))), 12L))
print(head(setdiff(rownames(exprs(eset)), raw[[1L]]), 12L))
