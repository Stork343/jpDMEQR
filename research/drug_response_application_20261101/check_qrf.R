library(jsonlite)
source("research/drug_response_application_20261101/comparator_prediction.R")
objects <- new.env()
load("results/scratch/empirical_selection_20261101/fbioselect_source/test_gdsc_data.RData", objects)
fd <- objects$fd
features <- list(y = fd$Response, genes = as.matrix(fd[, -(1:5)]), dose = fd$DosageSQ,
                   clinical = model.matrix(~ DRUG_NAME, fd), cluster_id = as.character(fd$ID))
first <- drug_qrf_fit(features, 1, 10L, 20261101L, trees = 100L)
second <- drug_qrf_fit(features, 1, 10L, 20261101L, trees = 100L)
predicted <- drug_qrf_predict(first, features, c(0.25, 0.5, 0.75))
repeated <- drug_qrf_predict(second, features, c(0.25, 0.5, 0.75))
stopifnot(all(is.finite(predicted)), max(abs(predicted - repeated)) == 0,
          all(predicted[, 1L] <= predicted[, 2L]), all(predicted[, 2L] <= predicted[, 3L]))
directory <- "results/preflight/empirical_selection_20261101/qrf"
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
write_json(list(status = "passed", actual_observations_used = TRUE,
                deterministic_error = max(abs(predicted - repeated)), quantiles_ordered = TRUE,
                bootstrap_unit = "experimental_id", package_version = as.character(packageVersion("ranger"))),
           file.path(directory, "checks.json"), auto_unbox = TRUE, pretty = TRUE)
cat("QRF 检查通过\n")
