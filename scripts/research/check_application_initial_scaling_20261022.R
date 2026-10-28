root <- normalizePath(".", mustWork = TRUE)
.libPaths(c(file.path(root, "results", "scratch", "profile_pair_package_20261022", "library"), .libPaths()))
library(jpDMEQR)
Sys.setenv(JPDMEQR_CLUSTER_SOLVER = "compiled")
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", environment())
source_v2_module(root, "application_v2", environment())
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(2L, 3L)) stop("Supply scaling, a new run ID and an optional saved source-design path.")
Sys.setenv(JPDMEQR_PROXIMAL_SCALING = args[1])
output <- file.path(root, "results", "scratch", args[2])
if (dir.exists(output)) stop("The predictive numerical validation directory already exists.")
dir.create(output, recursive = TRUE)
if (length(args) == 3L) {
  cached <- readRDS(args[3])
  metadata <- cached$metadata
  training <- cached$training
  design <- cached$design
} else {
cohort <- readRDS(file.path(root, "data", "derived", "GSE65391", "audit_20261022_bootstrap_v1", "analysis_cohort.rds"))
metadata <- cohort$eligible_metadata
folds <- cohort$subject_folds
heldout <- folds$subject_id[folds$repetition == 1L & folds$fold == 1L]
training <- !metadata$subject_id %in% heldout
rule <- yaml::read_yaml(file.path(root, "config", "module_scoring_20261022.yml"))
members <- read_frozen_module_score_members_v2(file.path(root, rule$score_membership_file))
signatures <- read_frozen_module_definitions_v2(file.path(root, rule$membership_file))
annotation <- read.csv(file.path(root, rule$annotation_path), stringsAsFactors = FALSE)
expression <- readRDS(cohort$expression_file)[, rownames(metadata), drop = FALSE]
design <- prepare_gse65391_prediction_design_v2(
  metadata, expression, annotation, signatures, members, rule, training,
  include_sledai = TRUE, filter_size = 5000L, screening_dimension = 100L, screening_tau = .5
)
}
saveRDS(list(design = design, metadata = metadata, training = training), file.path(output, "source_design.rds"))
started <- proc.time()[[3]]
fit <- withCallingHandlers(
  fit_application_profile_initial_v2(
    metadata$log_c3[training], design$X[training, , drop = FALSE],
    design$Z[training, , drop = FALSE], metadata$subject_id[training], tau = .5,
    penalty_factor = design$penalty_factor
  ),
  error = function(condition) {
    saveRDS(list(condition = condition, runtime_sec = proc.time()[[3]] - started),
            file.path(output, "failure_condition.rds"))
  }
)
saveRDS(fit, file.path(output, "initial_fit.rds"))
fit_runtime <- proc.time()[[3]] - started
reference <- fit_profile_lasso_v2(
  metadata$log_c3[training], design$X[training, , drop = FALSE],
  design$Z[training, , drop = FALSE], metadata$subject_id[training], .5,
  fit$h_est, fit$lambda_beta, lambda_gamma = fit$lambda_gamma,
  penalty_factor = fit$penalty_factor, beta_init = fit$beta,
  control = list(max_iter = 2000L, beta_tol = 1e-7, kkt_normalized_tol = 1e-3,
                  proximal_scaling = "none",
                  nuisance_control = list(grad_tol = 1e-8))
)
saveRDS(reference, file.path(output, "reference_warm_fixed_penalty.rds"))
if (!reference$converged || max(abs(reference$beta - fit$beta)) > 1e-6 ||
    abs(reference$composite_objective - fit$composite_objective) > 1e-8) {
  stop("The real predictive final fixed-penalty numerical comparison failed.")
}
write.csv(data.frame(
  p = ncol(design$X), n_subjects = length(unique(metadata$subject_id[training])),
  scaling = args[1], runtime_sec = fit_runtime,
  converged = fit$converged, kkt = fit$kkt_normalized,
  beta_change = fit$first_stage_beta_change, iterations = fit$iterations,
  max_reference_coefficient_difference = max(abs(reference$beta - fit$beta)),
  reference_objective_difference = abs(reference$composite_objective - fit$composite_objective),
  scope = "prediction_training_fold_numerical_computation_validation"
), file.path(output, "validation_summary.csv"), row.names = FALSE)
