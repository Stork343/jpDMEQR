source("research/submission_completion_20261027/load_comparators.R")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Provide a recorded SQR setting and replication.")
id <- args[1L]
replicate <- as.integer(args[2L])
path <- file.path("results/raw/simulation/submission_20261028_main_sqr/replications",
  id, sprintf("replicate_%05d/result.rds", replicate))
record <- readRDS(path)
config <- record$config
data <- prepare_analysis_data_v2(simulate_from_config_v2(config, record$seed), config, record$seed)$data
data$X <- cbind(`(Intercept)` = 1, data$X)
data$p <- ncol(data$X)
cv <- record$method_diagnostics[["SQR-DEBIASED-IID"]]$tuning$gamma_cv
fold <- 1L
training <- subset_clusters_v2(data, cv$roles[[fold]]$training_subjects)
validation <- subset_clusters_v2(data, cv$roles[[fold]]$validation_subjects)
initial <- cv$fold_initial_estimators[[fold]]
if (!identical(colnames(training$X), names(initial$beta)) ||
    nrow(training$X) != initial$N_training || nrow(validation$X) != initial$N_validation)
  stop("The reconstructed SQR training design differs.")
H <- sqr_iid_hessian_fn_v2(training$y, training$X, initial$beta, config$tau, initial$h)
validation_H <- sqr_iid_hessian_fn_v2(validation$y, validation$X, initial$beta, config$tau, initial$h)
cells <- cv$cells[cv$cells$fold == fold, ]
for (i in which(cells$status == "ok")) {
  saved <- cv$precision_records[[i]]$result
  coordinates <- match(rownames(saved$omega_raw_targets), colnames(training$X))
  canonical <- diag(ncol(training$X))[coordinates, , drop = FALSE]
  residual <- max(abs(saved$omega_raw_targets %*% H - canonical))
  if (!is.finite(residual) || residual > cells$gamma[i] + 1e-6)
    stop("The recorded CLIME directions do not match the reconstructed Hessian.")
}
output <- file.path("results/scratch/sqr_acceleration_20261030", paste0(id, "_", replicate))
dir.create(output, recursive = TRUE, showWarnings = FALSE)
saveRDS(list(H = H, validation_H = validation_H, X = training$X,
  training = training, validation = validation, initial = initial,
  cells = cells, precision_records = cv$precision_records[seq_len(nrow(cells))],
  reference_path = path, config = config, seed = record$seed),
  file.path(output, "inputs.rds"), compress = FALSE)
cat("Reconstructed", id, replicate, "matrix dimension", nrow(H), "\n")
