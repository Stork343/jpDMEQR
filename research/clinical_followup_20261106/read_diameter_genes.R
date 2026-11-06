library(GEOquery)
library(data.table)
src <- "data/raw/ispy_20261106"
pam <- fread(file.path(src, "pam50_model.csv"), skip = '"probe","probe.centroids"', colClasses = "character")
stopifnot(nrow(pam) == 50L, !anyDuplicated(pam$probe), !anyDuplicated(pam$EntrezGene.ID))
second <- fread(file.path(src, "expression.tsv.gz"), skip = 1L, header = FALSE, sep = "\t", select = 1L)[[1L]]
panel <- pam[probe %in% second]
platforms <- c("1708", "4133")
measures <- list()
rosters <- list()
for (id in platforms) {
  atlas <- fread(file.path(src, paste0("gpl", id, "_annot.gz")), skip = "ID\tGene title", header = TRUE, sep = "\t",
                 select = c("ID", "Gene ID"), fill = TRUE, colClasses = "character")
  array <- getGEO(filename = file.path(src, paste0("ispy1_gpl", id, ".gz")), getGPL = FALSE)
  subjects <- as.character(pData(array)[["i-spy id:ch2"]])
  keep <- !is.na(subjects) & grepl("^[0-9]+$", subjects)
  stopifnot(!anyDuplicated(subjects[keep]))
  rosters[[id]] <- data.table(patient = subjects[keep], platform = paste0("GPL", id))
  expr <- exprs(array)[, keep, drop = FALSE]
  colnames(expr) <- subjects[keep]
  atlas <- atlas[ID %in% rownames(expr) & `Gene ID` %in% pam$EntrezGene.ID]
  stopifnot(!anyDuplicated(atlas$ID))
  available <- pam[EntrezGene.ID %in% atlas[["Gene ID"]]]
  panel <- panel[probe %in% available$probe]
  out <- t(vapply(available$EntrezGene.ID, function(entrez) {
    selected <- atlas[`Gene ID` == entrez, ID]
    apply(expr[selected, , drop = FALSE], 2L, median, na.rm = TRUE)
  }, numeric(ncol(expr))))
  rownames(out) <- available$probe
  colnames(out) <- subjects[keep]
  measures[[id]] <- out
}
stopifnot(nrow(panel) >= 40L, all(c("MKI67", "BIRC5", "CCNB1") %in% panel$probe))
genes <- sort(panel$probe)
raw <- do.call(cbind, lapply(measures, function(x) x[genes, , drop = FALSE]))
stopifnot(!anyDuplicated(colnames(raw)))
fwrite(rbindlist(rosters), file.path(src, "ispy1_expression_roster.csv"))
fwrite(data.table(gene = genes, raw), file.path(src, "ispy1_panel_expression.csv"))
fwrite(panel[order(probe)], file.path(src, "common_pam50.csv"))
cat("Patients", ncol(raw), "shared genes", length(genes), "\n")
