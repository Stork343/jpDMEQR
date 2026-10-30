Sys.setenv(R_THREADS = "1")
source("research/score_validity_20261026/load_joint.R")
source_v2_module(".", "application_v2", .GlobalEnv)
source("research/submission_completion_20261027/dispatch_state.R")
source("research/joint_profile_manuscript_20261026/fixed_inference.R")
source("research/empirical_recovery_20261030/inputs.R")
input <- empirical_inputs()
processing <- "official_R1_R2_normexp_training_quantile"
raw <- readRDS("data/interim/GSE65391/raw_calibrated_20261030.rds")
if (!identical(colnames(raw$expression), rownames(input$metadata)))
  stop("The raw sensitivity visits differ from the registered cohort.")
training <- empirical_prediction_roles(input, 1L, 1L)$training[[1L]]
base <- empirical_prepare_base(input, training, processing)
empirical_validate_base(base, input$metadata, training, processing)
direct_target <- rowMeans(apply(raw$expression[, training, drop = FALSE], 2L, sort))
target <- base$processing_transformation$target
target_error <- max(abs(direct_target - target)) / max(1, max(abs(direct_target)))
if (target_error > 1e-12 || !setequal(base$processing_transformation$target_subjects,
    input$metadata$subject_id[training])) stop("The normalisation target differs from the training arrays.")
normalised <- preprocessCore::normalize.quantiles.use.target(raw$expression, target)
testing <- setdiff(seq_len(nrow(input$metadata)), training)[1L]
single <- preprocessCore::normalize.quantiles.use.target(raw$expression[, testing, drop = FALSE], target)
array_error <- max(abs(single[, 1L] - normalised[, testing])) / max(1, max(abs(single)))
if (array_error > 1e-12) stop("The heldout array transformation depends on other heldout arrays.")
full <- empirical_prepare_base(input, seq_len(nrow(input$metadata)), processing, filter_size = 0L)
X <- cbind(full$clinical$A$X, full$module_scores)
fit <- jp_fit_fixed(input$metadata$log_c3, X, matrix(1, nrow(X), 1L),
  input$metadata$subject_id, .5, coordinates = colnames(full$module_scores))
output <- "results/preflight/empirical_recovery_20261030"
submission_atomic_rds(list(fit = fit, target_error = target_error, array_error = array_error),
  "results/scratch/empirical_recovery_20261030/raw_check.rds")
jsonlite::write_json(list(passed = TRUE, processing = processing, subjects = 129L, visits = 867L,
  normalisation_target_relative_error = target_error, heldout_array_relative_error = array_error,
  full_sample_module_coordinates = nrow(fit$table), preprocessCore = as.character(packageVersion("preprocessCore")),
  checked = as.character(Sys.time())), file.path(output, paste0(processing, "_manifest.json")),
  auto_unbox = TRUE, pretty = TRUE)
cat("Verified the registered raw-array processing and full-sample module fit\n")
