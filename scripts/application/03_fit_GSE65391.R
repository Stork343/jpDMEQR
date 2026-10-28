root <- normalizePath(".", mustWork = TRUE)
.libPaths(c(file.path(root, "results", "scratch", "profile_pair_package_20261022", "library"), .libPaths()))
library(jpDMEQR)
Sys.setenv(JPDMEQR_CLUSTER_SOLVER = "compiled")
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", envir = environment())
source_v2_module(root, "application_v2", envir = environment())
source(file.path(root, "scripts", "application", "_paired_method_gate_20261022.R"))
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}
cohort_dir <- get_arg("cohort-dir", file.path(root, "data", "derived", "GSE65391", "audit_20261022_bootstrap_v1"))
cohort <- readRDS(file.path(cohort_dir, "analysis_cohort.rds"))
modules <- validate_gse65391_fit_inputs_v2(cohort$audit, file.path(root, "config", "module_signatures.csv"))
implementation <- validate_application_paired_gate_20261022(
  root, get_arg("checks-manifest", file.path(root, "results", "preflight", "application_checks.json")),
  method_version = "profile_pair_20261022"
)
rule_path <- file.path(root, "config", "module_scoring_20261022.yml")
rule <- yaml::read_yaml(rule_path)
members <- read_frozen_module_score_members_v2(file.path(root, rule$score_membership_file))
metadata <- cohort$eligible_metadata
validate_gse65391_subject_folds_v2(metadata, cohort$subject_folds)
repetition <- as.integer(get_arg("repetition", "1"))
fold <- as.integer(get_arg("fold", "1"))
model <- toupper(get_arg("model", "A"))
tau <- as.numeric(get_arg("tau", "0.5"))
filter_size <- as.integer(get_arg("features", "5000"))
screen_dimension <- as.integer(get_arg("screen-dimension", "100"))
lambda_multiplier <- as.numeric(get_arg("lambda-multiplier", "1"))
lambda_gamma <- as.numeric(get_arg("lambda-gamma", "1"))
h_multiplier <- as.numeric(get_arg("h-multiplier", "1"))
random_slope <- identical(get_arg("random-slope", "false"), "true")
if (!model %in% c("A", "B") || !filter_size %in% c(2000L, 5000L, 10000L) ||
    !screen_dimension %in% c(100L, 200L, 400L)) stop("Invalid frozen prediction design.")
folds <- cohort$subject_folds
test_subjects <- as.character(folds$subject_id[folds$repetition == repetition & folds$fold == fold])
if (!length(test_subjects)) stop("The requested outer fold is absent.")
train <- !as.character(metadata$subject_id) %in% test_subjects
expression <- readRDS(cohort$expression_file)[, rownames(metadata), drop = FALSE]
annotation <- read.csv(file.path(root, rule$annotation_path), stringsAsFactors = FALSE)
design <- prepare_gse65391_prediction_design_v2(
  metadata, expression, annotation, modules, members, rule, train,
  include_sledai = model == "B", filter_size = filter_size,
  screening_dimension = screen_dimension, screening_tau = tau, random_slope = random_slope
)
run_id <- get_arg("run-id", paste0("prediction_", format(Sys.time(), "%Y%m%dT%H%M%S")))
output <- file.path(root, "results", "raw", "application", run_id)
if (dir.exists(output)) stop("APPLICATION_OUTPUT_EXISTS: choose a new run_id.")
dir.create(output, recursive = TRUE)
saveRDS(list(clinical = design$clinical, module = design$gene_preprocessing,
             highdim = design$highdim_preprocessing, module_audit = design$module_audit),
        file.path(output, "preprocessing.rds"))
writeLines(implementation$commit_sha, file.path(output, "implementation_commit.txt"))
started <- proc.time()[[3]]
fit <- fit_application_profile_initial_v2(
  metadata$log_c3[train], design$X[train, , drop = FALSE], design$Z[train, , drop = FALSE],
  metadata$subject_id[train], tau,
  h_est = h_multiplier * (log(max(sum(design$penalty_factor > 0), 2L)) /
                            length(unique(metadata$subject_id[train])))^(1 / 4),
  lambda_gamma = lambda_gamma,
  lambda_multiplier = lambda_multiplier, penalty_factor = design$penalty_factor
)
runtime <- proc.time()[[3]] - started
saveRDS(fit, file.path(output, "fit.rds"))
prediction <- as.numeric(design$X[!train, , drop = FALSE] %*% fit$beta_l1)
outcome <- metadata$log_c3[!train]
loss <- (outcome - prediction) * (tau - as.numeric(outcome < prediction))
cluster_loss <- tapply(loss, metadata$subject_id[!train], mean)
cluster_calibration <- tapply(as.numeric(outcome <= prediction), metadata$subject_id[!train], mean)
train_sizes <- table(metadata$subject_id[train])
baseline <- quantreg::rq(metadata$log_c3[train] ~ 1, tau = tau,
                          weights = 1 / as.numeric(train_sizes[as.character(metadata$subject_id[train])]))
baseline_prediction <- as.numeric(stats::coef(baseline))
baseline_loss <- (outcome - baseline_prediction) * (tau - as.numeric(outcome < baseline_prediction))
baseline_mean <- mean(tapply(baseline_loss, metadata$subject_id[!train], mean))
metrics <- data.frame(
  run_id = run_id, repetition = repetition, fold = fold, model = model, tau = tau,
  method_version = fit$method_version, score_definition = fit$score_definition,
  filter_size = filter_size, screen_dimension = screen_dimension,
  lambda_multiplier = lambda_multiplier, lambda_gamma = lambda_gamma, h_multiplier = h_multiplier,
  h_est = fit$h_est, random_slope = random_slope,
  n_train_subjects = length(unique(metadata$subject_id[train])), n_train_visits = sum(train),
  n_test_subjects = length(test_subjects), n_test_visits = sum(!train),
  pinball_loss = mean(cluster_loss), intercept_pinball_loss = baseline_mean,
  intercept_skill = 1 - mean(cluster_loss) / baseline_mean,
  quantile_calibration = mean(cluster_calibration), selected_genes = sum(abs(fit$beta_l1[design$penalty_factor > 0]) > 1e-10),
  runtime_sec = runtime, status = "ok", converged = fit$converged,
  kkt_residual = fit$kkt_residual, profile_identity_error = fit$components$profile_identity_error,
  max_nuisance_gradient = fit$components$max_nuisance_gradient, stringsAsFactors = FALSE
)
write.csv(metrics, file.path(output, "fold_metrics.csv"), row.names = FALSE)
write.csv(data.frame(sample_id = rownames(metadata)[!train], subject_id = metadata$subject_id[!train],
                     observed_log_c3 = outcome, predicted_quantile = prediction, loss = loss),
          file.path(output, "predictions.csv"), row.names = FALSE)
write.csv(data.frame(coordinate = names(fit$beta_l1), beta_hat = as.numeric(fit$beta_l1),
                     selected = abs(fit$beta_l1) > 1e-10),
          file.path(output, "coefficients.csv"), row.names = FALSE)
jsonlite::write_json(list(run_id = run_id, implementation_commit = implementation$commit_sha,
                          numeric_certificate_commit = implementation$numeric_certificate_commit,
                          actual_source_commit = implementation$actual_source_commit,
                          solver_configuration = implementation$solver_configuration,
                          method_version = fit$method_version,
                          score_definition = fit$score_definition,
                          completed_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)),
                     file.path(output, "completed_manifest.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
print(metrics)
