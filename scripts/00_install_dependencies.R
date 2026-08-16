#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", hit[1])
}
as_bool <- function(x) tolower(as.character(x)) %in% c("true", "1", "yes", "y")

install_bioc <- as_bool(get_arg("bioc", "true"))
install_optional <- as_bool(get_arg("optional", "true"))

cran_required <- c(
  "quantreg", "CVXR", "testthat", "data.table", "yaml", "ggplot2",
  "dplyr", "tidyr", "readr", "purrr", "matrixStats", "digest", "lqmm"
)
cran_optional <- c("pbapply", "future.apply", "qs", "arrow", "peakRAM")

repos <- getOption("repos")
if (is.null(repos) || repos[["CRAN"]] == "@CRAN@") {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

install_missing_cran <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    message("Installing CRAN packages: ", paste(missing, collapse = ", "))
    install.packages(missing, dependencies = TRUE)
  }
  invisible(missing)
}

install_missing_cran(cran_required)
if (install_optional) install_missing_cran(cran_optional)

if (install_bioc) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  bioc <- c("GEOquery", "Biobase", "limma", "AnnotationDbi", "illuminaHumanv4.db")
  missing_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_bioc)) {
    message("Installing Bioconductor packages: ", paste(missing_bioc, collapse = ", "))
    BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
  }
}

message("Dependency check complete.")
print(sessionInfo())
