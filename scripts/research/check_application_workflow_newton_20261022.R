root <- normalizePath(".", mustWork = TRUE)
.libPaths(c(file.path(root, "results", "scratch", "profile_pair_package_20261022", "library"), .libPaths()))
library(jpDMEQR)
Sys.setenv(JPDMEQR_CLUSTER_SOLVER = "compiled",
           JPDMEQR_UNPENALISED_OPTIMIZER = "newton_armijo")
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", environment())
source_v2_module(root, "application_v2", environment())
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop("Supply model, q, error family, tau and a new numerical validation run ID.")
model <- args[1]
q <- as.integer(args[2])
error_family <- args[3]
tau <- as.numeric(args[4])
output <- file.path(root, "results", "scratch", args[5])
if (dir.exists(output)) stop("The complete-workflow numerical validation directory already exists.")
dir.create(output, recursive = TRUE)
config <- yaml::read_yaml(file.path(root, "config", "application_matched_prospective_20261022.yml"))
population <- readRDS(file.path(root, "data", "derived", "GSE65391", "execution_20261022_v1",
                                 sprintf("matched_population_model_%s_q%d.rds", model, q)))
rule <- yaml::read_yaml(file.path(root, "config", "module_scoring_20261022.yml"))
members <- read_frozen_module_score_members_v2(file.path(root, rule$score_membership_file))
signatures <- read_frozen_module_definitions_v2(file.path(root, rule$membership_file))
annotation <- read.csv(file.path(root, rule$annotation_path), stringsAsFactors = FALSE)
expression <- readRDS(file.path(root, "data", "derived", "GSE65391", "module_measurements_20261022_v1",
                                 "observed_module_probe_expression.rds"))
seed <- config$base_seed + (q + 2L * (match(error_family, config$error_families) - 1L)) * 1000000L + 1L
draw <- resample_application_matched_workflow_v2(population, expression, config, error_family, seed)
started <- proc.time()[[3]]
result <- withCallingHandlers(
  lowdim_profile_pair_workflow_inference_v2(
    draw$metadata, draw$expression, annotation, signatures, members, rule,
    tau = tau, h = 129^(-3 / 10), include_sledai = model == "B",
    random_slope = q == 2L, checkpoint_dir = file.path(output, "workflow_checkpoints"),
    checkpoint_storage = "compact"
  ),
  error = function(condition) saveRDS(condition, file.path(output, "failure_condition.rds"))
)
elapsed <- proc.time()[[3]] - started
diagnostics <- rbind(result$jackknife$full_sample$diagnostics, result$jackknife$diagnostics)
if (nrow(diagnostics) != 130L || any(!diagnostics$converged) ||
    any(diagnostics$refit_gradient > 1e-7) ||
    any(diagnostics$nuisance_gradient > 1e-7) ||
    any(diagnostics$n_coordinates != ncol(population$X)) ||
    any(result$table$jackknife_completed != 129L)) {
  stop("The declared complete 129-subject workflow numerical checks failed.")
}
save_compact_application_workflow_v2(result, file.path(output, "workflow_result.rds"))
replayed <- read_compact_application_workflow_v2(file.path(output, "workflow_result.rds"))
if (!isTRUE(all.equal(replayed$table, result$table, tolerance = 0, check.attributes = TRUE))) {
  stop("The complete workflow compact readback changed the inference table.")
}
write.csv(diagnostics, file.path(output, "all_fit_diagnostics.csv"), row.names = FALSE)
summary <- data.frame(
  model = model, q = q, error_family = error_family, tau = tau, seed = seed,
  n_subjects = 129L, n_visits = nrow(draw$metadata), p = ncol(population$X),
  completed_fits = nrow(diagnostics), completed_jackknives = 129L,
  runtime_sec = elapsed, max_refit_gradient = max(diagnostics$refit_gradient),
  max_nuisance_gradient = max(diagnostics$nuisance_gradient),
  max_identity_error = max(diagnostics$profile_identity_error),
  min_paired_curvature = min(diagnostics$hessian_min_eigenvalue),
  max_condition_number = max(diagnostics$hessian_condition_number),
  compact_readback_identical = TRUE,
  artifact_bytes = sum(file.info(list.files(output, recursive = TRUE, full.names = TRUE))$size),
  scope = "complete_workflow_numerical_computation_validation", stringsAsFactors = FALSE
)
write.csv(summary, file.path(output, "validation_summary.csv"), row.names = FALSE)
print(summary)
