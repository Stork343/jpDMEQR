#!/usr/bin/env Rscript
# Merge shard run directories (pilot_v3_*_s1..sK) into a single run directory.
args <- commandArgs(trailingOnly = TRUE)
root <- "D:/OneDrive/paper/Jointly Penalised and Debiased Inference for High-Dimensional Mixed-Effects Quantile Regression/jpDMEQR"
run_id <- args[1] %||% stop("merged run_id required")
n_shards <- if (length(args) > 1 && nzchar(args[2])) as.integer(args[2]) else stop("n_shards required")
raw_dir <- file.path(root, "results", "raw", "simulation")
merged <- file.path(raw_dir, run_id)
if (dir.exists(merged)) stop("Merged run directory already exists: ", merged)
dir.create(merged)

csvs <- c("replication_metrics.csv", "coordinate_metrics.csv", "theory_diagnostics.csv",
          "screening_records.csv", "runner_failures.csv")
for (cc in csvs) {
  parts <- lapply(seq_len(n_shards), function(s) {
    p <- file.path(raw_dir, paste0(run_id, "_s", s), cc)
    if (file.exists(p)) read.csv(p, stringsAsFactors = FALSE) else NULL
  })
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (length(parts)) {
    tab <- do.call(rbind, parts)
    if (nrow(tab)) {
      ord <- do.call(order, c(tab[intersect(c("experiment_id", "replicate", "method_id", "coordinate"),
                                             names(tab))], list(na.last = TRUE)))
      tab <- tab[ord, , drop = FALSE]
    }
    write.csv(tab, file.path(merged, cc), row.names = FALSE)
    cat(cc, ": rows =", nrow(tab), "\n")
  }
}
# provenance
meta <- lapply(seq_len(n_shards), function(s) {
  p <- file.path(raw_dir, paste0(run_id, "_s", s))
  list(shard = s, run_dir = p,
       commit = if (file.exists(file.path(p, "implementation_commit.txt")))
         readLines(file.path(p, "implementation_commit.txt"), warn = FALSE)[1] else NA_character_,
       config_sha = if (file.exists(file.path(p, "config_sha256.txt")))
         readLines(file.path(p, "config_sha256.txt"), warn = FALSE)[1] else NA_character_)
})
writeLines(jsonlite::toJSON(meta, auto_unbox = TRUE, pretty = TRUE),
           file.path(merged, "merge_manifest.json"))
cat("MERGED:", merged, "\n")
