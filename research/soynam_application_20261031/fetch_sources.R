library(data.table)
library(jsonlite)
sources <- data.table(package = c("NAM", "SoyNAM"),
  version = c("1.8.0", "1.6.2"),
  revision = c("e2980d8fc7be457e0a5e65abead7b96a1ce01111", "34b1c01912ba6b95b67b78d98667b50d0602904b"),
  object_path = c("data/met.RData", "data/soyin.RData"),
  destination = c("met.RData", "soyin_reference.RData"))
destination <- "data/raw/soynam_20261031"
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
transport <- Sys.getenv("SOYNAM_GIT_TRANSPORT", "https")
if (!transport %in% c("https", "ssh")) stop("SOYNAM_GIT_TRANSPORT must be https or ssh.")
for (i in seq_len(nrow(sources))) {
  output <- file.path(destination, sources$destination[i])
  if (file.exists(output)) next
  checkout <- file.path("results/scratch/soynam_20261031", paste0(sources$package[i], "_source"))
  url <- if (transport == "ssh") paste0("ssh://git@ssh.github.com:443/cran/", sources$package[i], ".git") else
    paste0("https://github.com/cran/", sources$package[i], ".git")
  if (!dir.exists(checkout)) {
    status <- system2("git", c("clone", "--no-checkout", shQuote(url), shQuote(checkout)))
    if (status != 0L) stop("The source repository download failed.")
  }
  status <- system2("git", c("-C", shQuote(checkout), "show",
    shQuote(paste0(sources$revision[i], ":", sources$object_path[i]))), stdout = output)
  if (status != 0L) stop("The pinned source object could not be extracted.")
}
met <- new.env()
reference <- new.env()
load(file.path(destination, "met.RData"), envir = met)
load(file.path(destination, "soyin_reference.RData"), envir = reference)
if (!all(c("Gen", "Obs") %in% ls(met)) || !all(c("data.line.in", "data.check.in") %in% ls(reference)))
  stop("The source object names differ from their registered release.")
source_sizes <- file.info(file.path(destination, sources$destination))$size
if (anyNA(source_sizes) || any(source_sizes <= 0)) stop("The source files are missing or empty.")
sources[, bytes := source_sizes]
fwrite(sources, "results/preflight/soynam_20261031/source_manifest.csv")
print(sources)
