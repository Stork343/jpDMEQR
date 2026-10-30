args <- commandArgs(trailingOnly = TRUE)
release <- args[1L]
if (length(args) != 1L || !release %in% c("R1", "R2")) stop("Supply R1 or R2.")
name <- paste0("GSE65391_non-normalized_data_Illumina_HT12_V4_", release, ".txt.gz")
url <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65391/suppl/", name)
directory <- "data/raw/GSE65391"
path <- file.path(directory, name)
manifest <- file.path(directory, paste0("download_recovery_20261030_", release, ".csv"))
if (file.exists(path)) stop("The requested raw expression file already exists.")
partial <- paste0(path, ".partial")
status <- system2("curl", c("--fail", "--location", "--connect-timeout", "30",
  "--max-time", "7200", "--continue-at", "-", "--output", shQuote(partial), shQuote(url)))
if (status != 0L) stop("The raw expression download did not complete.")
if (system2("gzip", c("-t", shQuote(partial))) != 0L) stop("The downloaded gzip file is incomplete.")
if (!file.rename(partial, path)) stop("Unable to publish the downloaded input.")
write.csv(data.frame(release = release, url = url, path = path,
  downloaded_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  bytes = file.info(path)$size, md5 = unname(tools::md5sum(path)), status = "downloaded"),
  manifest, row.names = FALSE)
cat("Downloaded and verified", release, file.info(path)$size, "bytes\n")
