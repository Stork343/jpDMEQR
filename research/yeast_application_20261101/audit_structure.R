library(data.table)
library(jsonlite)
source <- new.env()
load("data/raw/yeast_bloom_2015/extracted/genotypes.RData", envir = source)
load("data/raw/yeast_bloom_2015/extracted/phenotypes.RData", envir = source)
print(ls(source))
if (!all(c("gdata", "pheno_raw") %in% ls(source))) stop("The expected genotype and replicate objects are absent.")
G <- source$gdata
phenotypes <- source$pheno_raw
cat("Genotype dimensions:", dim(G), "\n")
cat("Genotype class:", class(G), "\n")
cat("First strain identifiers:\n")
print(head(rownames(G), 24L))
cat("Last strain identifiers:\n")
print(tail(rownames(G), 12L))
cat("First marker identifiers:\n")
print(head(colnames(G), 12L))
cat("Phenotype names:\n")
print(names(phenotypes))
summary <- rbindlist(lapply(names(phenotypes), function(trait) {
  current <- phenotypes[[trait]]
  data.table(trait, object_class = paste(class(current), collapse = "/"),
    strains = length(current), identifiers = length(names(current)),
    observations = sum(lengths(current)),
    minimum_replicates = min(lengths(current)), maximum_replicates = max(lengths(current)),
    finite_observations = sum(vapply(current, function(x) sum(is.finite(x)), integer(1))),
    complete_pairs = sum(vapply(current, function(x) length(x) == 2L && all(is.finite(x)), logical(1))),
    missing_genotype_identifiers = sum(!names(current) %in% rownames(G)))
}))
destination <- "results/preflight/empirical_selection_20261101/yeast_bloom_2015"
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
fwrite(summary, file.path(destination, "measurement_structure.csv"))
write_json(list(source_url = "https://static-content.springer.com/esm/art%3A10.1038%2Fncomms9712/MediaObjects/41467_2015_BFncomms9712_MOESM729_ESM.zip",
  paper_doi = "10.1038/ncomms9712", genotype_dimensions = dim(G),
  phenotype_count = length(phenotypes), genotype_missing = sum(is.na(G)),
  strain_identifier_examples = head(rownames(G), 24L),
  marker_identifier_examples = head(colnames(G), 12L),
  outcome_associations_inspected = FALSE), file.path(destination, "structure.json"),
  auto_unbox = TRUE, pretty = TRUE)
print(summary)
