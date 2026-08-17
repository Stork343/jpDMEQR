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
  "quantreg",
  "CVXR",
  "clarabel",
  "ECOSolveR",
  "scs",
  "testthat",
  "data.table",
  "yaml",
  "jsonlite",
  "digest",
  "ggplot2",
  "dplyr",
  "tidyr",
  "readr",
  "purrr",
  "matrixStats",
  "lqmm"
)
cran_optional <- c(
  "pbapply",
  "future.apply",
  "parallelly",
  "qs",
  "arrow",
  "peakRAM"
)

repos <- getOption("repos")
if (is.null(repos) || is.na(repos[["CRAN"]]) || repos[["CRAN"]] == "@CRAN@") {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

install_missing_cran <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    message("Installing CRAN packages: ", paste(missing, collapse = ", "))
    install.packages(missing, dependencies = TRUE)
  }
  still_missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still_missing)) {
    stop(
      "Required CRAN packages remain unavailable: ",
      paste(still_missing, collapse = ", ")
    )
  }
  invisible(missing)
}

install_missing_cran(cran_required)
if (install_optional) {
  optional_missing <- cran_optional[
    !vapply(cran_optional, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(optional_missing)) {
    message(
      "Installing optional CRAN packages: ",
      paste(optional_missing, collapse = ", ")
    )
    try(
      install.packages(optional_missing, dependencies = TRUE),
      silent = FALSE
    )
  }
}

if (install_bioc) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  bioc_required <- c(
    "GEOquery",
    "Biobase",
    "limma",
    "AnnotationDbi",
    "illuminaHumanv4.db"
  )
  missing_bioc <- bioc_required[
    !vapply(bioc_required, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_bioc)) {
    message(
      "Installing Bioconductor packages: ",
      paste(missing_bioc, collapse = ", ")
    )
    BiocManager::install(
      missing_bioc,
      ask = FALSE,
      update = FALSE
    )
  }
  still_missing_bioc <- bioc_required[
    !vapply(bioc_required, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(still_missing_bioc)) {
    stop(
      "Required Bioconductor packages remain unavailable: ",
      paste(still_missing_bioc, collapse = ", ")
    )
  }
}

solvers <- CVXR::installed_solvers()
approved <- intersect(c("CLARABEL", "ECOS", "SCS"), solvers)
if (!length(approved)) {
  stop(
    "CVXR is installed but no approved solver is available. Installed solvers: ",
    paste(solvers, collapse = ", ")
  )
}

message("Dependency check complete.")
message("Approved CVXR solvers: ", paste(approved, collapse = ", "))
print(sessionInfo())
