args <- commandArgs(trailingOnly = TRUE)
library(filelock)
release <- args[1L]
if (length(args) != 1L || !release %in% c("R1", "R2")) stop("Supply R1 or R2.")
name <- paste0("GSE65391_non-normalized_data_Illumina_HT12_V4_", release, ".txt.gz")
url <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65391/suppl/", name)
directory <- "data/raw/GSE65391"
path <- file.path(directory, name)
lock <- filelock::lock(file.path(directory, paste0(release, "_acquisition.lock")), timeout = 0)
if (is.null(lock)) stop("This raw release is already being downloaded.")
manifest <- file.path(directory, paste0("download_recovery_20261030_", release, ".csv"))
if (file.exists(path)) stop("The requested raw expression file already exists.")
partial <- paste0(path, ".pypdl.partial")
logs <- file.path("results/scratch/empirical_recovery_20261030", paste0("download_", release))
dir.create(logs, recursive = TRUE, showWarnings = FALSE)
environment <- Sys.getenv()
environment[["PYTHONPATH"]] <- normalizePath("results/scratch/python_dependencies")
processx::run("python3", c(normalizePath("research/empirical_recovery_20261030/download_file.py"),
  url, file.path(normalizePath(dirname(partial)), basename(partial))), wd = normalizePath(logs), env = environment,
  stdout = file.path(normalizePath(logs), "transfer.log"), stderr = "2>&1")
if (system2("gzip", c("-t", shQuote(partial))) != 0L) stop("The downloaded gzip file is incomplete.")
if (!file.rename(partial, path)) stop("Unable to publish the downloaded input.")
write.csv(data.frame(release = release, url = url, path = path,
  downloaded_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  bytes = file.info(path)$size, md5 = unname(tools::md5sum(path)), status = "downloaded"),
  manifest, row.names = FALSE)
cat("Downloaded and verified", release, file.info(path)$size, "bytes\n")
filelock::unlock(lock)
