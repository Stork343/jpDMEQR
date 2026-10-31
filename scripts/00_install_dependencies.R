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
  "quantreg", "CVXR", "clarabel", "testthat", "pkgload", "withr",
  "data.table", "yaml", "ps",
  "jsonlite", "ggplot2", "dplyr", "tidyr", "readr", "readxl", "purrr",
  "matrixStats", "digest", "lqmm", "Matrix", "quadprog", "nloptr",
  "Rcpp", "RcppArmadillo", "RcppEigen", "glmnet", "filelock", "processx", "rcmdcheck", "statmod", "geeVerse",
  "osqp", "igraph", "knitr", "kableExtra", "hqreg", "mgcv"
)
cran_optional <- c(
  "pbapply", "future.apply", "qs", "arrow", "peakRAM", "foreach",
  "doParallel", "ECOSolveR", "ncvreg", "osqp"
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
  remaining <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(remaining)) stop("Required CRAN packages remain unavailable: ", paste(remaining, collapse = ", "))
  invisible(missing)
}

install_missing_cran(cran_required)
install_highs_interface <- function() {
  destination <- file.path("results", "scratch", "dependencies", "library")
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  installed <- utils::installed.packages(lib.loc = destination)
  if (!"highs" %in% rownames(installed) ||
      package_version(installed["highs", "Version"]) < "1.15.1-1") {
    old_timeout <- getOption("timeout")
    options(timeout = max(600, old_timeout))
    on.exit(options(timeout = old_timeout), add = TRUE)
    install.packages("highs", lib = destination, type = "source")
  }
  if (utils::packageVersion("highs", lib.loc = destination) < "1.15.1-1")
    stop("The required HiGHS basis interface is unavailable.")
  .libPaths(c(normalizePath(destination, mustWork = TRUE), .libPaths()))
}
install_highs_interface()
install_parametric_clime <- function() {
  if (requireNamespace("fastclime", quietly = TRUE)) {
    if (packageVersion("fastclime") != "1.4.1.1") stop("The parametric CLIME interface requires fastclime 1.4.1.1.")
    return(invisible(NULL))
  }
  scratch <- file.path("results", "scratch", "dependencies")
  dir.create(scratch, recursive = TRUE, showWarnings = FALSE)
  scratch <- normalizePath(scratch, mustWork = TRUE)
  archive <- file.path(scratch, "fastclime_1.4.1.1.tar.gz")
  old_timeout <- getOption("timeout")
  options(timeout = max(600, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)
  download.file("https://cran.r-project.org/src/contrib/Archive/fastclime/fastclime_1.4.1.1.tar.gz",
    archive, mode = "wb")
  install.packages(archive, repos = NULL, type = "source")
  if (!requireNamespace("fastclime", quietly = TRUE) || packageVersion("fastclime") != "1.4.1.1")
    stop("The required archived CLIME dependency was not installed.")
}
install_parametric_clime()
source("research/submission_completion_20261027/install_compiled_clime.R")
if (install_optional) {
  optional_missing <- cran_optional[!vapply(cran_optional, requireNamespace, logical(1), quietly = TRUE)]
  if (length(optional_missing)) {
    message("Installing optional CRAN packages: ", paste(optional_missing, collapse = ", "))
    try(install.packages(optional_missing, dependencies = TRUE), silent = TRUE)
  }
}

if (install_bioc) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  bioc <- c("GEOquery", "Biobase", "limma", "preprocessCore", "AnnotationDbi", "illuminaHumanv4.db")
  missing_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_bioc)) {
    message("Installing Bioconductor packages: ", paste(missing_bioc, collapse = ", "))
    BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
  }
}

solvers <- CVXR::installed_solvers()
approved <- intersect(c("CLARABEL", "ECOS", "SCS"), solvers)
if (!length(approved)) {
  stop("CVXR is installed but no approved solver (CLARABEL/ECOS/SCS) is available.")
}
message("Approved CVXR solvers: ", paste(approved, collapse = ", "))
message("Dependency check complete.")
print(sessionInfo())
