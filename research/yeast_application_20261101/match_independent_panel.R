library(data.table)
library(tidyr)
library(jsonlite)

original <- fread("data/raw/yeast_bloom_2013/BYxRM_PhenoData.txt")
stopifnot(nrow(original) == 1008L, !anyDuplicated(original$V1))
expression <- new.env()
growth <- new.env()
load("results/scratch/empirical_selection_20261101/eqtl_source/RData/counts.RData", expression)
load("data/raw/yeast_bloom_2015/extracted/genotypes.RData", growth)
first <- expression$counts$gdata
second <- growth$gdata
stopifnot(!anyDuplicated(rownames(first)), !anyDuplicated(rownames(second)))
first_markers <- as.data.table(separate_wider_delim(
  data.frame(marker = colnames(first)), marker, delim = "_", names = c("coordinate", "alleles")))
first_markers <- as.data.table(separate_wider_delim(first_markers, coordinate,
  delim = ":", names = c("chromosome", "position")))
first_markers <- as.data.table(separate_wider_delim(first_markers, alleles,
  delim = "/", names = c("BY", "RM")))
first_markers[, first_index := seq_len(.N)]
second_markers <- as.data.table(separate_wider_delim(data.frame(marker = colnames(second)), marker,
  delim = "_", names = c("genome_position", "chromosome", "position", "BY", "RM")))
second_markers[, second_index := seq_len(.N)]
shared <- merge(first_markers, second_markers, by = c("chromosome", "position", "BY", "RM"))
stopifnot(nrow(shared) >= 10000L, uniqueN(shared$chromosome) == 16L)
shared[, position := as.integer(position)]
setorder(shared, chromosome, position)
anchor <- shared[, .SD[unique(as.integer(round(seq(1, .N, length.out = 20L))))], by = chromosome]
matched <- data.table(panel_row = seq_len(nrow(original)), strain_id = original$V1)
matched[, expression_2018_available := strain_id %in% rownames(first)]
available <- which(matched$expression_2018_available)
similarity <- tcrossprod(first[matched$strain_id[available], anchor$first_index, drop = FALSE],
                         second[, anchor$second_index, drop = FALSE]) / nrow(anchor)
best <- max.col(similarity, ties.method = "first")
second_best <- apply(similarity, 1L, function(values) sort(values, decreasing = TRUE)[2L])
matched[, `:=`(growth_2015_id = NA_character_, anchor_concordance = NA_real_,
                second_anchor_concordance = NA_real_)]
matched[available, `:=`(growth_2015_id = rownames(second)[best],
  anchor_concordance = (similarity[cbind(seq_along(best), best)] + 1) / 2,
  second_anchor_concordance = (second_best + 1) / 2)]
old_genotypes <- first[matched$strain_id[available], shared$first_index, drop = FALSE]
new_genotypes <- second[matched$growth_2015_id[available], shared$second_index, drop = FALSE]
matched[, complete_marker_concordance := NA_real_]
matched[available, complete_marker_concordance := rowMeans(old_genotypes == new_genotypes)]
matched[, verified_2015 := !is.na(complete_marker_concordance) & complete_marker_concordance >= 0.995 &
                          anchor_concordance - second_anchor_concordance >= 0.05]
stopifnot(!anyDuplicated(matched[verified_2015 == TRUE, growth_2015_id]))
output_dir <- "results/preflight/empirical_selection_20261101/yeast_panel_match"
derived_dir <- "data/derived/yeast_20261101"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(matched, file.path(derived_dir, "independent_panel_match.csv"))
fwrite(matched[, .(rows = .N), by = .(expression_2018_available, verified_2015)],
       file.path(output_dir, "match_quality.csv"))
write_json(list(status = "genotype_matching_completed", original_panel_rows = 1008L,
                 shared_markers = nrow(shared), verified_2015_rows = sum(matched$verified_2015),
                 minimum_verified_concordance = min(matched[verified_2015 == TRUE, complete_marker_concordance]),
                 outcome_associations_used = FALSE,
                 independent_panel_source = "10.1038/nature11867",
                 expression_source = "10.7554/eLife.35471"),
           file.path(output_dir, "checks.json"), auto_unbox = TRUE, pretty = TRUE)
print(matched[, .N, by = .(expression_2018_available, verified_2015)])
print(summary(matched$complete_marker_concordance))
