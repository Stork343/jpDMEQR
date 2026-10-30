library(ps)
library(filelock)
source("research/empirical_recovery_20261030/execution_control.R")
empirical_require_active()
root <- "results/scratch/empirical_recovery_20261030"
lock <- filelock::lock(file.path(root, "raw_preparation.lock"), timeout = 0)
if (is.null(lock)) stop("The raw preparation is already running.")
downloads <- read.csv(file.path(root, "active_downloads.csv"))
if (!setequal(downloads$release, c("R1", "R2"))) stop("Both registered raw releases must have an acquisition process.")
repeat {
  downloads <- read.csv(file.path(root, "active_downloads.csv"))
  finished <- vapply(downloads$release, function(release) {
    file.exists(file.path("data/raw/GSE65391", paste0("download_recovery_20261030_", release, ".csv")))
  }, logical(1))
  if (all(finished)) break
  if (any(!finished & !downloads$pid %in% ps::ps_pids())) stop("A raw-array acquisition ended before producing its complete source record.")
  Sys.sleep(10)
}
Sys.setenv(R_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1", OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
if (!file.exists("data/interim/GSE65391/raw_calibrated_20261030.rds")) {
  status <- system2(file.path(R.home("bin"), "Rscript"), "research/empirical_recovery_20261030/prepare_raw.R")
  if (status != 0L) stop("The raw-array preparation did not pass.")
}
status <- system2(file.path(R.home("bin"), "Rscript"), "research/empirical_recovery_20261030/check_raw.R")
if (status != 0L) stop("The raw-array processing checks did not pass.")
destination <- file.path(path.expand("~"), "jpDMEQR-backups", "data/raw/GSE65391")
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
for (release in c("R1", "R2")) {
  filename <- paste0("GSE65391_non-normalized_data_Illumina_HT12_V4_", release, ".txt.gz")
  if (!file.copy(file.path("data/raw/GSE65391", filename), file.path(destination, filename), overwrite = TRUE, copy.date = TRUE))
    stop("The verified source backup was not written.")
  manifest <- paste0("download_recovery_20261030_", release, ".csv")
  if (!file.copy(file.path("data/raw/GSE65391", manifest), file.path(destination, manifest), overwrite = TRUE, copy.date = TRUE))
    stop("The acquisition record backup was not written.")
}
cat("The registered raw-array workflows are ready for the empirical dispatcher\n")
filelock::unlock(lock)
