library(data.table)
library(yaml)
config <- read_yaml("config/application_milk_20261031.yml")
source("research/score_validity_20261026/load_joint.R")
source("research/milk_application_20261031/nodewise.R")
cohorts <- readRDS(file.path(config$results_directory, "cohorts.rds"))
cache_path <- file.path(config$results_directory, "design_cache.rds")
if (file.exists(cache_path)) stop("The verified design cache already exists.")
cache <- list()
checks <- list()
for (cohort in names(cohorts)) {
  current <- cohorts[[cohort]]
  data <- sb_data(current$y, current$X, current$Z, current$cluster_id)
  assignment <- pc_folds(data$rows, 2L)
  for (fold in 1:2) {
    training <- sb_subset(data, unlist(data$rows[assignment != fold], use.names = FALSE))
    factor <- training$X * sqrt(training$weight)
    started <- proc.time()[[3L]]
    fit <- milk_nodewise_compute(factor, training$rows, 1L)
    seconds <- proc.time()[[3L]] - started
    label <- paste(cohort, fold, sep = "_")
    cache[[label]] <- list(factor = factor, rows = training$rows,
      subjects = names(training$rows), coordinate = 1L, fit = fit)
    checks[[label]] <- data.table(cohort = cohort, fold = fold, clusters = training$n,
      path_kkt = max(fit$tuning$path_kkt), selected_kkt = fit$kkt,
      normalisation_error = fit$normalisation_error, identity_error = fit$identity_error,
      selected_ratio = fit$tuning$selected_ratio, runtime_seconds = seconds)
    print(checks[[label]])
  }
}
entry <- cache[[1L]]
columns <- c(1:13, seq(14L, ncol(entry$factor), length.out = 30L))
columns <- unique(as.integer(columns))
reference <- pc_nodewise(entry$factor[, columns, drop = FALSE], entry$rows, 1L)
checked <- milk_nodewise_compute(entry$factor[, columns, drop = FALSE], entry$rows, 1L)
coefficient_error <- max(abs(reference$coefficient - checked$coefficient))
loss_error <- max(abs(reference$tuning$fold_loss_sum - checked$tuning$fold_loss_sum))
if (coefficient_error > 1e-7 || loss_error > 1e-7 ||
    reference$tuning$selected != checked$tuning$selected) stop("The actual-design reference check failed.")
saveRDS(cache, cache_path, compress = FALSE)
fwrite(rbindlist(checks), "results/preflight/milk_scs_20261031/design_cache_checks.csv")
jsonlite::write_json(list(coefficient_error = coefficient_error, validation_loss_error = loss_error,
  matching_selected_ratio = TRUE, actual_observations = nrow(entry$factor),
  reference_columns = length(columns)), "results/preflight/milk_scs_20261031/nodewise_reference_check.json",
  auto_unbox = TRUE, pretty = TRUE, digits = NA)
if (!file.copy(cache_path, config$backup_directory, overwrite = FALSE)) stop("Design cache backup failed.")
cat("Verified actual-design conditional regressions and saved their shared cache.\n")
