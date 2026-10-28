required <- c("Rcpp", "RcppArmadillo", "RcppEigen", "glmnet", "filelock", "processx")
available <- vapply(required, requireNamespace, logical(1), quietly = TRUE)
if (!all(available)) stop(paste("Install the declared dependencies:", paste(required[!available], collapse = ", ")))
library_path <- "results/scratch/score_validity_20261026/library"
dir.create(library_path, recursive = TRUE, showWarnings = FALSE)
status <- system2(file.path(R.home("bin"), "R"),
                   c("CMD", "INSTALL", "--no-multiarch", paste0("--library=", shQuote(library_path)), "."))
if (status != 0L) stop("Local package installation failed.")
description <- utils::packageDescription("jpDMEQR", lib.loc = library_path)
saveRDS(list(library = normalizePath(library_path), package = description, session = sessionInfo()),
         "results/scratch/score_validity_20261026/installation.rds")
cat("Installed jpDMEQR in", normalizePath(library_path), "\n")
