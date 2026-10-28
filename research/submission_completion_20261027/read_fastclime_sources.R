directory <- "results/scratch/submission_sources_20261027"
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
archive <- file.path(directory, "fastclime_1.4.1.1.tar.gz")
download.file("https://cran.r-project.org/src/contrib/Archive/fastclime/fastclime_1.4.1.1.tar.gz", archive, mode = "wb")
untar(archive, exdir = directory)
cat(archive, "\n")
