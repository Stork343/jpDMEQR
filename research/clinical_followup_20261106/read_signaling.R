library(data.table)
library(GEOquery)
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 1L, args %in% c("roster", "discovery"))
path <- "data/raw/ispy_20261106"
study <- getGEO(filename=file.path(path,"protein_matrix.txt.gz"), GSEMatrix=TRUE, getGPL=FALSE, AnnotGPL=FALSE)
patients <- as.character(pData(study)$title)
stopifnot(length(patients) == 736L, !anyDuplicated(patients), nrow(exprs(study)) == 139L)
if (args == "roster") {
  roster <- data.table(patient=patients)
  fwrite(roster, file.path(path,"signaling_roster.csv"))
} else {
  inputs <- "data/derived/ispy2_signaling_20261106"
  roles <- fread(file.path(inputs,"roles.csv"), colClasses=list(character="patient"))
  ids <- roles[role=="discovery",patient]
  selected <- match(ids, patients)
  stopifnot(!anyNA(selected))
  signals <- t(exprs(study)[,selected,drop=FALSE])
  stopifnot(!anyDuplicated(colnames(signals)), all(is.finite(signals[!is.na(signals)])))
  fwrite(cbind(data.table(patient=ids),as.data.table(signals)), file.path(inputs,"discovery_signals.csv"))
}
