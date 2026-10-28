root <- normalizePath(".", mustWork = TRUE)
.libPaths(c(file.path(root, "results", "scratch", "profile_pair_package_20261022", "library"), .libPaths()))
library(jpDMEQR)
Sys.setenv(JPDMEQR_CLUSTER_SOLVER = "compiled")
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", envir = environment())
source_v2_module(root, "application_v2", envir = environment())
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}
config <- yaml::read_yaml(file.path(root, "config", "application_matched_prospective_20261022.yml"))
implementation <- validate_application_population_implementation_v2(root, get_arg("core-build-commit", "5c47a62"))
source_identity <- validate_application_numeric_source_identity_v2(
  root, get_arg("application-build-commit", Sys.getenv("JPDMEQR_APPLICATION_BUILD_COMMIT", unset = ""))
)
registry <- read.csv(file.path(root, "data", "derived", "GSE65391", "execution_20261022_v1", "matched_registry.csv"),
                      stringsAsFactors = FALSE)
case <- registry[registry$case_id == get_arg("case-id"), , drop = FALSE]
repeat_id <- as.integer(get_arg("replicate", "1"))
asset_path <- get_arg("target-asset")
if (nrow(case) != 1L || repeat_id < 1L || repeat_id > config$replications || is.null(asset_path)) {
  stop("Supply a complete declared matched case, replicate and target asset.")
}
asset <- readRDS(asset_path)
expected_case <- sprintf("AM_%s_q%d_%s", case$model, case$q, case$error_family)
if (!identical(asset$case, expected_case) || !identical(asset$method_version, "profile_pair_20261022") ||
    !identical(asset$config, config) ||
    length(asset$target_records) != 4L || any(!vapply(asset$target_records, function(record) {
      record$n_population >= 100000L && isTRUE(record$lower$converged) && isTRUE(record$refinement$converged)
    }, logical(1)))) stop("The complete matching population target asset is required.")
population <- readRDS(file.path(root, "data", "derived", "GSE65391", "execution_20261022_v1",
                                 sprintf("matched_population_model_%s_q%d.rds", case$model, case$q)))
if (!application_population_reference_equal_v2(
  asset$population_reference, population[c("X", "Z", "rows", "beta0", "preprocessing")]
)) {
  stop("The target population design, measurement definitions or working parameterisation changed.")
}
rule <- yaml::read_yaml(file.path(root, "config", "module_scoring_20261022.yml"))
members <- read_frozen_module_score_members_v2(file.path(root, rule$score_membership_file))
signatures <- read_frozen_module_definitions_v2(file.path(root, rule$membership_file))
annotation <- read.csv(file.path(root, rule$annotation_path), stringsAsFactors = FALSE)
expression <- readRDS(file.path(root, "data", "derived", "GSE65391", "module_measurements_20261022_v1",
                                 "observed_module_probe_expression.rds"))
seed <- config$base_seed + (case$q + 2L * (match(case$error_family, config$error_families) - 1L)) * 1000000L + repeat_id
draw <- resample_application_matched_workflow_v2(population, expression, config, case$error_family, seed)
run_id <- get_arg("run-id", "matched_20261022")
output <- file.path(root, "results", "raw", "application_matched", run_id, case$case_id,
                     sprintf("replicate_%04d", repeat_id))
if (dir.exists(output)) stop("The matched replicate directory already exists.")
dir.create(output, recursive = TRUE)
saveRDS(list(case = case, repeat_id = repeat_id, seed = seed, target_asset = asset_path,
             source_subject_types = draw$generated$source_subject_types,
             implementation_commit = implementation$commit_sha,
             source_identity = source_identity,
             source_inputs = list(
               population = file.path("data", "derived", "GSE65391", "execution_20261022_v1",
                                       sprintf("matched_population_model_%s_q%d.rds", case$model, case$q)),
               expression = file.path("data", "derived", "GSE65391", "module_measurements_20261022_v1",
                                       "observed_module_probe_expression.rds"),
               annotation = rule$annotation_path, module_members = rule$score_membership_file,
               module_signatures = rule$membership_file,
               module_rule = "config/module_scoring_20261022.yml",
               generation_config = "config/application_matched_prospective_20261022.yml"
             ),
             solver_configuration = solver_configuration_v2()), file.path(output, "run_request.rds"))
record_failure <- function(condition) {
    target <- asset$targets[asset$targets$tau == case$tau, ]
    failure <- data.frame(case_id = case$case_id, replicate = repeat_id,
                            coordinate = target$coordinate, target = target$target,
                            beta_hat = NA_real_, sandwich_se = NA_real_, jackknife_se = NA_real_,
                            status = "failed", failure_stage = condition$failure_stage %||% "workflow_postprocessing",
                            failure_message = conditionMessage(condition), stringsAsFactors = FALSE)
    write.csv(failure, file.path(output, "coordinate_metrics_failed.csv"), row.names = FALSE)
    saveRDS(condition, file.path(output, "failure_condition.rds"))
}
withCallingHandlers({
result <- lowdim_profile_pair_workflow_inference_v2(
  draw$metadata, draw$expression, annotation, signatures, members, rule, tau = case$tau,
  h = config$n_subjects^(-3 / 10), include_sledai = case$model == "B",
  random_slope = case$q == 2L, checkpoint_dir = file.path(output, "workflow_checkpoints"),
  checkpoint_storage = "compact"
)
target <- asset$targets[asset$targets$tau == case$tau, ]
mapped <- application_matched_coordinate_metrics_v2(result, population, target, case, repeat_id)
workflow_path <- file.path(output, "workflow_result.rds")
save_compact_application_workflow_v2(
  result, workflow_path,
  extra = list(population_unit_map = mapped$population_unit_map,
               mapped_jackknife = mapped$mapped_jackknife, target_asset = asset_path,
               case = case, repeat_id = repeat_id)
)
replayed <- read_compact_application_workflow_v2(workflow_path)
replayed_metrics <- application_matched_coordinate_metrics_v2(replayed, population, target, case, repeat_id)
if (!isTRUE(all.equal(replayed_metrics$table, mapped$table, tolerance = 0, check.attributes = TRUE))) {
  stop("The compact matched workflow did not reproduce the complete population-unit reporting table.")
}
write.csv(mapped$table, file.path(output, "coordinate_metrics.csv"), row.names = FALSE)
jsonlite::write_json(list(case_id = case$case_id, replicate = repeat_id, seed = seed,
                          implementation_commit = implementation$commit_sha,
                          source_identity = source_identity,
                          core_build_commit = implementation$core_build_commit,
                          population_target_build_commit = asset$implementation_commit,
                          method_version = "profile_pair_20261022", completed_subject_jackknives = 129L,
                          solver_configuration = solver_configuration_v2(),
                          storage_mode = "complete_workflow_compact_v1",
                          target_asset = normalizePath(asset_path),
                          completed_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)),
                     file.path(output, "completed_manifest.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
}, error = record_failure)
