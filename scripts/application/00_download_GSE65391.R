#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
`%||%` <- function(x, y) if (is.null(x)) y else x
parse_args <- function(args) {
  out <- list()
  for (x in args) {
    if (!grepl("^--", x)) next
    z <- strsplit(sub("^--", "", x), "=", fixed = TRUE)[[1]]
    out[[z[1]]] <- if (length(z) > 1) paste(z[-1], collapse = "=") else "true"
  }
  out
}
as_bool <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  tolower(as.character(x)) %in% c("true", "1", "yes", "y")
}
cli <- parse_args(args)

full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "scripts/application/00_download_GSE65391.R"
root <- if (file.exists("config/application.yml")) "." else normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)

processed <- as_bool(cli$processed, TRUE)
raw_tar <- as_bool(cli$raw, FALSE)
non_normalized <- as_bool(cli$`non-normalized`, FALSE)
force <- as_bool(cli$force, FALSE)
retries <- as.integer(cli$retries %||% 3L)
out_dir <- cli$`out-dir` %||% file.path(root, "data", "raw", "GSE65391")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

urls <- c(
  series_matrix = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65391/matrix/GSE65391_series_matrix.txt.gz",
  raw_tar = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65391/suppl/GSE65391_RAW.tar",
  non_normalized_r1 = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65391/suppl/GSE65391_non-normalized_data_Illumina_HT12_V4_R1.txt.gz",
  non_normalized_r2 = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65391/suppl/GSE65391_non-normalized_data_Illumina_HT12_V4_R2.txt.gz"
)
selected <- character()
if (processed) selected <- c(selected, "series_matrix")
if (raw_tar) selected <- c(selected, "raw_tar")
if (non_normalized) selected <- c(selected, "non_normalized_r1", "non_normalized_r2")
if (!length(selected)) stop("Nothing selected for download.")

url_basename <- function(url) basename(sub("[?#].*$", "", url))
download_one <- function(name, url) {
  dest <- file.path(out_dir, url_basename(url))
  if (file.exists(dest) && file.info(dest)$size > 0 && !force) {
    message("Using existing file: ", dest)
    return(data.frame(name = name, url = url, local_path = normalizePath(dest),
                      downloaded_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
                      bytes = file.info(dest)$size, md5 = unname(tools::md5sum(dest)),
                      status = "existing", message = "", stringsAsFactors = FALSE))
  }
  last_error <- ""
  for (attempt in seq_len(retries)) {
    tmp <- paste0(dest, ".part")
    if (file.exists(tmp)) unlink(tmp)
    ok <- tryCatch({
      utils::download.file(url, tmp, mode = "wb", quiet = FALSE,
                           method = if (capabilities("libcurl")) "libcurl" else "auto")
      file.exists(tmp) && file.info(tmp)$size > 0
    }, error = function(e) {
      last_error <<- conditionMessage(e)
      FALSE
    })
    if (isTRUE(ok)) {
      if (!file.rename(tmp, dest)) stop("Cannot move completed download to ", dest)
      return(data.frame(name = name, url = url, local_path = normalizePath(dest),
                        downloaded_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
                        bytes = file.info(dest)$size, md5 = unname(tools::md5sum(dest)),
                        status = "downloaded", message = "", stringsAsFactors = FALSE))
    }
    Sys.sleep(min(30, 2^attempt))
  }
  data.frame(name = name, url = url, local_path = dest,
             downloaded_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
             bytes = NA_real_, md5 = NA_character_, status = "failed",
             message = last_error, stringsAsFactors = FALSE)
}

manifest <- do.call(rbind, lapply(selected, function(nm) download_one(nm, urls[[nm]])))
manifest_path <- file.path(out_dir, "download_manifest.csv")
write.csv(manifest, manifest_path, row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out_dir, "download_sessionInfo.txt"))
print(manifest)
if (any(manifest$status == "failed")) stop("One or more downloads failed; inspect ", manifest_path)
message("Download manifest: ", manifest_path)
