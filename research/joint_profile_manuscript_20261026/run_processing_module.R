source("research/score_validity_20261026/load_joint.R")
source("research/joint_profile_manuscript_20261026/fixed_inference.R")
source_v2_module(".", "application_v2", .GlobalEnv)
source("research/submission_completion_20261027/dispatch_state.R")
source("research/empirical_recovery_20261030/inputs.R")
args <- commandArgs(trailingOnly = TRUE)
processing <- args[1L]
model <- args[2L]
tau <- as.numeric(args[3L])
if (!processing %in% c("official_R1_R2_normexp_training_quantile", "processed_median_across_mapped_probes") ||
    !model %in% c("A", "B") || !tau %in% c(.25, .5, .75)) stop("Invalid processing module specification.")
input <- empirical_inputs()
metadata <- input$metadata
subjects <- unique(as.character(metadata$subject_id))
n <- length(subjects)
cache <- file.path("data/derived/GSE65391/empirical_recovery_20261030/preprocessing", processing)
output <- file.path("results/raw/application/joint_profile_processing_20261026", processing,
                     paste0(model, "_tau", round(100 * tau)))
dir.create(output, recursive = TRUE, showWarnings = FALSE)
lock <- filelock::lock(file.path(output, "run.lock"), timeout = 0)
if (is.null(lock)) stop("This processing module workflow is already running.")
if (file.exists(file.path(output, "result.rds"))) quit(status = 0L)
request_path <- file.path(output, "request.rds")
if (file.exists(request_path)) {
  prior <- readRDS(request_path)
  if (prior$processing != processing || prior$model != model || prior$tau != tau || prior$n_subjects != n)
    stop("The module workflow definition changed.")
} else submission_atomic_rds(list(processing = processing, model = model, tau = tau, source_cache = cache,
  n_subjects = n, implementation_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  source_inputs = input$source, session = sessionInfo(), started = Sys.time()), request_path)
globalCallingHandlers(error = function(condition) saveRDS(condition, file.path(output, "failure.rds")))
full <- NULL
estimates <- matrix(NA_real_, n, 5L)
omitted <- character(n)
diagnostics <- list()
for (i in 0:n) {
  base <- empirical_processing_module_base(input, processing,
    omitted = if (i == 0L) NA_character_ else subjects[i], index = i)
  position <- match(base$sample_ids, rownames(metadata))
  if (anyNA(position) || anyDuplicated(position)) stop("Processing sample identifiers differ.")
  rows <- base$training_rows
  md <- metadata[position[rows], , drop = FALSE]
  included <- unique(as.character(md$subject_id))
  if (i == 0L) {
    if (!setequal(included, subjects)) stop("Incomplete full processing cohort.")
  } else {
    missing <- setdiff(subjects, included)
    if (length(missing) != 1L) stop("A deletion stage must omit exactly one subject.")
    omitted[i] <- missing
  }
  module_scores <- base$module_scores
  colnames(module_scores) <- paste0("module__", sub("^module__", "", colnames(module_scores)))
  X <- cbind(base$clinical[[model]]$X, module_scores)[rows, , drop = FALSE]
  targets <- colnames(module_scores)
  checkpoint <- file.path(output, "fits", sprintf("fit_%03d.rds", i + 1L))
  if (file.exists(checkpoint)) {
    stored <- readRDS(checkpoint)
    if (!identical(stored$sample_ids, rownames(md)) || stored$model != model || stored$tau != tau ||
        stored$processing != processing) stop("A saved module fit has a different definition.")
    fitted <- stored$fit
  } else {
    fitted <- jp_fit_fixed(md$log_c3, X, matrix(1, nrow(X), 1L), md$subject_id,
      tau, coordinates = targets, h = n^(-.3))
    submission_atomic_rds(list(fit = fitted, sample_ids = rownames(md), model = model,
      tau = tau, processing = processing), checkpoint)
  }
  if (i == 0L) {
    full <- fitted
    colnames(estimates) <- targets
  } else {
    if (!identical(targets, colnames(estimates))) stop("Processing module coordinates changed.")
    estimates[i, ] <- fitted$beta[targets]
  }
  diagnostics[[i + 1L]] <- fitted$diagnostics
}
if (!setequal(omitted, subjects) || anyDuplicated(omitted) || any(!is.finite(estimates))) stop("Incomplete processing workflow.")
centered <- sweep(estimates, 2L, colMeans(estimates), "-")
table <- full$table
table$workflow_se <- sqrt((n - 1) / n * colSums(centered^2))
table$processing <- processing
table$model <- model
table$tau <- tau
table$family_size <- 30L
critical <- qnorm(1 - .05 / (2 * table$family_size))
table$simultaneous_lower <- table$beta_hat - critical * table$se
table$simultaneous_upper <- table$beta_hat + critical * table$se
table$workflow_simultaneous_lower <- table$beta_hat - critical * table$workflow_se
table$workflow_simultaneous_upper <- table$beta_hat + critical * table$workflow_se
submission_atomic_rds(list(table = table, full = full, estimates = estimates, omitted_subjects = omitted,
             diagnostics = do.call(rbind, diagnostics), completed = Sys.time()), file.path(output, "result.rds"))
write.csv(table, file.path(output, "estimates.csv"), row.names = FALSE)
filelock::unlock(lock)
