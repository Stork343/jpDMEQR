library(data.table)
library(Biobase)
library(limma)
library(filelock)
source("research/submission_completion_20261027/dispatch_state.R")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L || (length(args) && !args %in% c("R1", "R2"))) stop("Supply one declared release or no arguments.")
releases <- if (length(args)) args else c("R1", "R2")
output <- "data/interim/GSE65391/raw_calibrated_20261030.rds"
if (file.exists(output)) stop("The calibrated raw expression already exists.")
cohort <- readRDS("data/derived/GSE65391/audit_20261022_bootstrap_v1/analysis_cohort.rds")
metadata <- cohort$eligible_metadata
eset <- readRDS("data/interim/GSE65391/gse65391_eset.rds")
phenotype <- Biobase::pData(eset)
observed_probes <- rownames(Biobase::exprs(eset))
if (!all(startsWith(phenotype$description, "Sample name: ")))
  stop("The official phenotype does not contain complete array identifiers.")
arrays <- data.frame(array_id = trimws(sub("^Sample name: ", "", phenotype$description)),
  sample_id = rownames(phenotype), subject_id = phenotype[["subject:ch1"]],
  disease = phenotype[["disease state:ch1"]], stringsAsFactors = FALSE)
if (anyDuplicated(arrays$array_id) || anyDuplicated(arrays$sample_id) || nrow(arrays) != 996L)
  stop("The official source array identifiers are incomplete.")
rm(eset, phenotype)
gc()
corrected <- checks <- list()
for (release in releases) {
  path <- file.path("data/raw/GSE65391", paste0("GSE65391_non-normalized_data_Illumina_HT12_V4_", release, ".txt.gz"))
  manifest <- read.csv(file.path("data/raw/GSE65391", paste0("download_recovery_20261030_", release, ".csv")))
  if (nrow(manifest) != 1L || manifest$bytes != file.info(path)$size || manifest$status != "downloaded")
    stop("The raw-array download does not match its completed acquisition record.")
  saved <- file.path("data/interim/GSE65391", paste0("raw_normexp_20261030_", release, ".rds"))
  if (file.exists(saved)) {
    value <- readRDS(saved)
    corrected[[release]] <- value$expression
    checks[[release]] <- value$checks
    next
  }
  raw <- fread(path, check.names = FALSE, nThread = 1L, showProgress = FALSE)
  signal_columns <- seq.int(2L, ncol(raw), by = 2L)
  detection_columns <- signal_columns + 1L
  if (names(raw)[1L] != "ID_REF" || ncol(raw) != 1L + 2L * length(signal_columns) ||
      any(names(raw)[detection_columns] != "Detection Pval") ||
      anyDuplicated(raw[[1L]]) || !setequal(raw[[1L]], observed_probes) ||
      anyNA(match(names(raw)[signal_columns], arrays$array_id)))
    stop("The raw-array expression and detection columns do not match the official records.")
  signal <- as.matrix(raw[, ..signal_columns])
  detection <- as.matrix(raw[, ..detection_columns])
  rownames(signal) <- rownames(detection) <- raw[[1L]]
  colnames(detection) <- colnames(signal)
  rm(raw)
  gc()
  if (any(!is.finite(signal)) || any(!is.finite(detection)) || any(detection < 0 | detection > 1))
    stop("The raw intensities or detection probabilities are incomplete.")
  fitted <- limma::nec(signal, detection.p = detection, offset = 16, robust = FALSE)
  single <- limma::nec(signal[, 1L, drop = FALSE], detection.p = detection[, 1L, drop = FALSE],
    offset = 16, robust = FALSE)
  difference <- max(abs(fitted[, 1L] - single[, 1L]))
  if (any(!is.finite(fitted)) || any(fitted <= 0) || difference > 1e-8)
    stop("The per-array background correction failed its independent check.")
  fitted <- fitted[match(observed_probes, rownames(fitted)), , drop = FALSE]
  checks[[release]] <- data.frame(release = release, probes = nrow(fitted), arrays = ncol(fitted),
    single_array_difference = difference)
  submission_atomic_rds(list(expression = fitted, checks = checks[[release]],
    source_file = path, correction = list(method = "limma::nec", offset = 16, robust = FALSE,
      package_version = as.character(packageVersion("limma")))), saved, compress = TRUE)
  corrected[[release]] <- fitted
  rm(signal, detection, fitted, single)
  gc()
}
if (length(releases) == 1L) {
  print(checks[[releases]])
  quit(status = 0L)
}
raw_ids <- unlist(lapply(corrected, colnames), use.names = FALSE)
if (anyDuplicated(raw_ids) || !setequal(raw_ids, arrays$array_id))
  stop("The two releases do not reconcile all observed arrays.")
arrays$release <- ifelse(arrays$array_id %in% colnames(corrected$R1), "R1", "R2")
healthy <- arrays[arrays$disease == "Healthy", ]
paired <- intersect(healthy$subject_id[healthy$release == "R1"], healthy$subject_id[healthy$release == "R2"])
if (length(paired) != 24L || length(intersect(paired, metadata$subject_id)))
  stop("The healthy calibration subjects do not match the declared technical pairs.")
pair_ids <- lapply(c("R1", "R2"), function(release) {
  rows <- healthy[healthy$release == release & healthy$subject_id %in% paired, ]
  if (nrow(rows) != length(paired) || anyDuplicated(rows$subject_id))
    stop("Every calibration subject must occur once in each release.")
  rows$array_id[match(paired, rows$subject_id)]
})
ratio <- rowMeans(corrected$R2[, pair_ids[[2L]], drop = FALSE]) /
  rowMeans(corrected$R1[, pair_ids[[1L]], drop = FALSE])
if (any(!is.finite(ratio)) || any(ratio <= 0)) stop("The healthy calibration ratios are invalid.")
corrected$R1 <- corrected$R1 * ratio
expression <- cbind(corrected$R1, corrected$R2)
colnames(expression) <- arrays$sample_id[match(colnames(expression), arrays$array_id)]
if (anyNA(match(rownames(metadata), colnames(expression)))) stop("A retained clinical visit is absent from the raw arrays.")
expression <- expression[, rownames(metadata), drop = FALSE]
submission_atomic_rds(list(expression = expression, healthy_ratio = ratio,
  calibration_subjects = paired, array_mapping = arrays, checks = do.call(rbind, checks),
  source = "official_GSE65391_R1_R2", correction = "limma_nec_offset16_robust_false",
  calibration = "per_probe_mean_R2_over_mean_R1_healthy_pairs", created = Sys.time()), output, compress = TRUE)
directory <- "results/preflight/empirical_recovery_20261030"
write.csv(do.call(rbind, checks), file.path(directory, "raw_array_checks.csv"), row.names = FALSE)
jsonlite::write_json(list(source = output, source_arrays = nrow(arrays),
  eligible_visits = ncol(expression), probes = nrow(expression), calibration_pairs = length(paired),
  ratio_range = range(ratio), limma_version = as.character(packageVersion("limma")),
  prepared = as.character(Sys.time())), file.path(directory, "raw_preparation_manifest.json"),
  auto_unbox = TRUE, pretty = TRUE)
cat("Prepared", ncol(expression), "observed raw-array visits with", length(paired), "healthy calibration pairs\n")
