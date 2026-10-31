library(data.table)
library(readxl)
library(jsonlite)

metadata <- as.data.table(read_excel("data/raw/yeast_bloom_2013/supplementary_table_1.xls"))
published <- metadata[Trait == "Manganese_Sulfate"]
stopifnot(nrow(published) == 1L, published$Dose == "10 mM")
objects <- new.env()
load("data/raw/yeast_bloom_2013/pheno_raw.RData", objects)
values <- objects$pheno_raw$Magnesium_Sulfate
counts <- vapply(values, function(value) sum(is.finite(value)), integer(1))
stopifnot(sum(counts >= 1L) == published[["Number of Segregants with at least one measurement"]],
          sum(counts >= 2L) == published[["Number of Segregants with at least 2 replicate measurements"]])
directory <- "results/preflight/empirical_selection_20261101/yeast_trait_identity"
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
mapping <- data.table(source_field = "Magnesium_Sulfate", published_trait = published$Trait,
                      chemical = "manganese sulfate", concentration = published$Dose,
                      measured_strains = sum(counts >= 1L), replicated_strains = sum(counts >= 2L),
                      reference = "10.1038/nature11867 Supplementary Table 1")
fwrite(mapping, file.path(directory, "trait_mapping.csv"))
write_json(list(status = "verified_from_source_table", source_label_preserved = TRUE,
                source_table = "data/raw/yeast_bloom_2013/supplementary_table_1.xls"),
           file.path(directory, "checks.json"), auto_unbox = TRUE, pretty = TRUE)
print(mapping)
