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
config <- yaml::read_yaml(file.path(root, "config", "application_execution_20261022.yml"))
implementation <- validate_application_paired_gate_20261022(
  root, get_arg("checks-manifest", file.path(root, "results", "preflight", "application_checks.json")),
  method_version = "profile_pair_20261022"
)
r <- as.integer(get_arg("repetition", "1"))
f <- as.integer(get_arg("fold", "1"))
model <- get_arg("model", "A")
tau <- as.numeric(get_arg("tau", "0.5"))
filter_size <- as.integer(get_arg("features", "5000"))
random_slope <- identical(get_arg("random-slope", "false"), "true")
method <- get_arg("method", "PROFILE-INITIAL")
if (!method %in% c("PROFILE-INITIAL", "POOLED-QR-LASSO", "QGEE-SCAD", "LQMM")) stop("Unknown prediction method.")
if (random_slope && method != "PROFILE-INITIAL") stop("The declared random-slope validation is for the profile method.")
if (!r %in% 1:5 || !f %in% 1:5 || !model %in% config$models ||
    !tau %in% config$quantiles || !filter_size %in% config$prediction$filter_sizes) stop("Invalid validation specification.")
obj <- readRDS(file.path(root, config$cohort_directory, "analysis_cohort.rds"))
modules <- validate_gse65391_fit_inputs_v2(obj$audit, file.path(root, "config", "module_signatures.csv"))
metadata <- obj$eligible_metadata
outer <- obj$subject_folds
heldout <- outer$subject_id[outer$repetition == r & outer$fold == f]
outer_training <- !metadata$subject_id %in% heldout
inner <- read.csv(file.path(root, "data", "derived", "GSE65391", "execution_20261022_v1",
                            "inner_subject_folds.csv"), stringsAsFactors = FALSE)
inner <- inner[inner$outer_repetition == r & inner$outer_fold == f, ]
if (length(intersect(inner$subject_id, heldout)) ||
    !setequal(inner$subject_id, metadata$subject_id[outer_training])) stop("The inner fold source is inconsistent.")
rule <- yaml::read_yaml(file.path(root, "config", "module_scoring_20261022.yml"))
members <- read_frozen_module_score_members_v2(file.path(root, rule$score_membership_file))
annotation <- read.csv(file.path(root, rule$annotation_path), stringsAsFactors = FALSE)
expression <- readRDS(obj$expression_file)[, rownames(metadata), drop = FALSE]
grid <- expand.grid(h_multiplier = config$prediction$h_multiplier_grid,
                     lambda_gamma = config$prediction$lambda_gamma_grid,
                     lambda_multiplier = config$prediction$lambda_multiplier_grid,
                     screen_dimension = config$prediction$screening_dimensions,
                     stringsAsFactors = FALSE)
if (method == "POOLED-QR-LASSO") {
  grid <- expand.grid(h_multiplier = 1, lambda_gamma = 1,
                       lambda_multiplier = config$comparators$pooled$lambda_multipliers,
                       screen_dimension = config$comparators$pooled$screening_dimensions)
} else if (method == "QGEE-SCAD") {
  grid <- expand.grid(h_multiplier = 1, lambda_gamma = 1, lambda_multiplier = 1,
                       screen_dimension = config$comparators$qgee$screening_dimensions)
} else if (method == "LQMM") {
  grid <- data.frame(h_multiplier = 1, lambda_gamma = 1, lambda_multiplier = 1, screen_dimension = 0L)
}
grid$candidate_id <- seq_len(nrow(grid))
run_id <- get_arg("run-id", sprintf("cv_%s_%s_r%d_f%d_tau%02d_filter%d", method, model, r, f, round(100 * tau), filter_size))
output <- file.path(root, "results", "raw", "application_validation", run_id)
resume <- identical(get_arg("resume", "false"), "true")
request <- list(numeric_certificate_commit = implementation$numeric_certificate_commit,
                 solver_configuration = implementation$solver_configuration,
                 initial_acceleration = Sys.getenv("JPDMEQR_INITIAL_ACCELERATION", unset = "none"),
                 r = r, f = f, model = model,
                 random_slope = random_slope,
                 method = method,
                 tau = tau, filter_size = filter_size, grid = grid, inner = inner)
if (dir.exists(output)) {
  if (!resume || file.exists(file.path(output, "completed_manifest.json")) ||
      !identical(readRDS(file.path(output, "run_request.rds")), request)) stop("Validation resume identity does not match.")
} else {
  dir.create(output, recursive = TRUE)
  saveRDS(request, file.path(output, "run_request.rds"))
}
withCallingHandlers({
scores <- list()
for (inner_fold in seq_len(5L)) {
  clinical_cache <- new.env(parent = emptyenv())
  validation_subjects <- inner$subject_id[inner$inner_fold == inner_fold]
  train <- outer_training & !metadata$subject_id %in% validation_subjects
  test <- metadata$subject_id %in% validation_subjects
  for (dimension in unique(grid$screen_dimension)) {
    design_path <- file.path(output, sprintf("preprocessing_inner%d_screen%d.rds", inner_fold, dimension))
    design <- if (file.exists(design_path)) readRDS(design_path) else {
      value <- if (method == "LQMM") {
        item <- prepare_gse65391_lowdim_design_v2(
          metadata, expression, annotation, modules, train, rule,
          include_sledai = model == "B", module_score_members = members
        )
        item$penalty_factor <- rep(0, ncol(item$X))
        item
      } else prepare_gse65391_prediction_design_v2(
          metadata, expression, annotation, modules, members, rule, train,
          include_sledai = model == "B", filter_size = filter_size,
          screening_dimension = dimension, screening_tau = tau, random_slope = random_slope
        )
      environment(value$clinical$formula) <- baseenv()
      value$clinical$X <- NULL
      saveRDS(value, design_path)
      value
    }
    for (index in which(grid$screen_dimension == dimension)) {
      candidate <- grid[index, ]
      path <- file.path(output, sprintf("inner%d_candidate%03d.rds", inner_fold, candidate$candidate_id))
      result <- if (file.exists(path)) readRDS(path) else {
        value <- fit_application_prediction_candidate_v2(
          metadata, design, train, test, tau, candidate$h_multiplier,
          candidate$lambda_gamma, candidate$lambda_multiplier, method = method,
          seed = config$base_seed + r * 100000L + f * 10000L + inner_fold * 100L,
          initial_control = list(clinical_cache = clinical_cache)
        )
        value$candidate <- candidate
        value$record$cluster_scores <- NULL
        saveRDS(value, path)
        value
      }
      result$subject_metrics$candidate_id <- candidate$candidate_id
      result$subject_metrics$inner_fold <- inner_fold
      scores[[length(scores) + 1L]] <- result$subject_metrics
    }
  }
}
scores <- do.call(rbind, scores)
cv <- aggregate(loss ~ candidate_id, scores, mean)
selected_id <- cv$candidate_id[order(cv$loss, cv$candidate_id)][1]
selected <- grid[grid$candidate_id == selected_id, ]
write.csv(scores, file.path(output, "inner_subject_losses.csv"), row.names = FALSE)
write.csv(merge(grid, cv), file.path(output, "inner_cv_candidates.csv"), row.names = FALSE)
design <- if (method == "LQMM") {
  item <- prepare_gse65391_lowdim_design_v2(
    metadata, expression, annotation, modules, outer_training, rule,
    include_sledai = model == "B", module_score_members = members
  )
  item$penalty_factor <- rep(0, ncol(item$X))
  item
} else prepare_gse65391_prediction_design_v2(
  metadata, expression, annotation, modules, members, rule, outer_training,
  include_sledai = model == "B", filter_size = filter_size,
  screening_dimension = selected$screen_dimension, screening_tau = tau, random_slope = random_slope
)
final <- fit_application_prediction_candidate_v2(
  metadata, design, outer_training, !outer_training, tau,
  selected$h_multiplier, selected$lambda_gamma, selected$lambda_multiplier, method = method,
  seed = config$base_seed + r * 100000L + f * 10000L
)
for (name in c("repetition", "fold", "model", "tau", "filter_size")) {
  value <- switch(name, repetition = r, fold = f, model = model, tau = tau, filter_size = filter_size)
  final$metrics[[name]] <- value
  final$predictions[[name]] <- value
}
final$metrics$run_id <- run_id
final$metrics$random_slope <- random_slope
saveRDS(final$record, file.path(output, "outer_fit.rds"))
saveRDS(design, file.path(output, "outer_preprocessing.rds"))
write.csv(final$metrics, file.path(output, "fold_metrics.csv"), row.names = FALSE)
write.csv(final$predictions, file.path(output, "predictions.csv"), row.names = FALSE)
write.csv(selected, file.path(output, "selected_tuning.csv"), row.names = FALSE)
jsonlite::write_json(list(run_id = run_id, implementation_commit = implementation$commit_sha,
                          numeric_certificate_commit = implementation$numeric_certificate_commit,
                          actual_source_commit = implementation$actual_source_commit,
                          solver_configuration = implementation$solver_configuration,
                          initial_acceleration = request$initial_acceleration,
                          method_version = final$record$method_version, method = method,
                          inner_folds = 5L, candidates = nrow(grid), selected_candidate = selected_id,
                          completed_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)),
                     file.path(output, "completed_manifest.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
print(final$metrics)
}, error = function(condition) {
  saveRDS(condition, file.path(output, "failure_condition.rds"))
  jsonlite::write_json(list(status = "failed", failure_stage = condition$failure_stage %||% "validation_workflow",
                            message = conditionMessage(condition), implementation = implementation),
                       file.path(output, "failed_manifest.json"), auto_unbox = TRUE, pretty = TRUE)
})
