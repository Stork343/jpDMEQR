library(data.table)
rows <- readRDS("data/derived/drugmatrix_20261101/discovery.rds")[exploration_role == "training"]
panel <- readRDS("data/derived/drugmatrix_20261101/panel.rds")
records <- list()
for (fold in 0:3) {
  current <- if (fold == 0L) rows else rows[inner_fold != fold]
  C <- panel$C[current$sample_id, , drop = FALSE]
  decomposition <- qr(C)
  columns <- sort(decomposition$pivot[seq_len(decomposition$rank)])
  for (id in unique(current$compound_id)) {
    retained_rank <- qr(C[current$compound_id != id, columns, drop = FALSE])$rank
    if (retained_rank < length(columns)) records[[length(records) + 1L]] <- data.table(
      training_fold = fold, compound = id, full_rank = length(columns), retained_rank = retained_rank,
      animals = sum(current$compound_id == id))
  }
}
output <- "results/preflight/empirical_selection_20261101/drugmatrix_offset_checks"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
fwrite(rbindlist(records), file.path(output, "background_design_ranks.csv"))
print(rbindlist(records))
