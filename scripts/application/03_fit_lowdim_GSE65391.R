full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else
  "scripts/application/03_fit_lowdim_GSE65391.R"
root <- if (file.exists("R/application_v2.R")) "." else
  normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
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
cohort_dir <- get_arg("cohort-dir")
if (is.null(cohort_dir)) stop("Supply --cohort-dir=<versioned audit directory>.")
obj <- readRDS(file.path(cohort_dir, "analysis_cohort.rds"))
modules <- validate_gse65391_fit_inputs_v2(obj$audit, file.path(root, "config", "module_signatures.csv"))
implementation <- validate_application_paired_gate_20261022(
  root, get_arg("checks-manifest", file.path(root, "results", "preflight", "application_checks.json")),
  inference = TRUE, method_version = "profile_pair_20261022"
)
module_rule_path <- get_arg("module-rule", file.path(root, "config", "module_scoring_20261022.yml"))
module_rule <- yaml::read_yaml(module_rule_path)
score_members_path <- get_arg("score-members", file.path(root, "config", "module_score_members_20261022.csv"))
score_members <- read_frozen_module_score_members_v2(score_members_path)
score_source <- "bandwidth_pair"
if (!identical(implementation$method_version, "profile_pair_20261022")) {
  stop("APPLICATION_METHOD_GATE_FAILED: a matching paired-method manifest is required.")
}
models <- strsplit(toupper(get_arg("models", "A,B")), ",", fixed = TRUE)[[1]]
if (!all(models %in% c("A", "B"))) stop("Invalid clinical adjustment model.")
taus <- as.numeric(strsplit(get_arg("quantiles", "0.25,0.50,0.75"), ",", fixed = TRUE)[[1]])
h_multiplier <- as.numeric(get_arg("h-multiplier", "1"))
lambda_gamma <- as.numeric(get_arg("lambda-gamma", "1"))
random_slope <- identical(get_arg("random-slope", "false"), "true")
if (!all(taus %in% c(0.10, 0.25, 0.50, 0.75)) ||
    !h_multiplier %in% c(0.75, 1, 1.25) || !lambda_gamma %in% c(0.5, 1, 2)) {
  stop("The requested low-dimensional specification is outside the frozen sensitivity design.")
}
if (0.10 %in% taus) {
  diagnostics <- obj$audit$curvature$diagnostics
  selected <- diagnostics$tau == 0.10 & abs(diagnostics$h_multiplier - h_multiplier) < 1e-12
  if (sum(selected) != 1L || !isTRUE(diagnostics$pass[selected])) {
    stop("CURVATURE_GATE_FAILED: the requested tau .10 bandwidth did not pass its audit.")
  }
}
run_id <- get_arg("run-id", paste0("gse65391_lowdim_", format(Sys.time(), "%Y%m%dT%H%M%S")))
out_dir <- file.path(root, "results", "raw", "application", run_id)
if (dir.exists(out_dir)) stop("APPLICATION_OUTPUT_EXISTS: choose a new run_id.")
metadata <- obj$eligible_metadata
expression <- readRDS(obj$expression_file)
expression_probe <- expression[, rownames(metadata), drop = FALSE]
annotation_path <- get_arg("annotation-path", file.path(root, module_rule$annotation_path))
annotation <- read.csv(annotation_path,
                       stringsAsFactors = FALSE, check.names = FALSE)
dir.create(out_dir, recursive = TRUE)
writeLines(implementation$commit_sha, file.path(out_dir, "implementation_commit.txt"))
file.copy(module_rule_path, file.path(out_dir, "module_rule.yml"), overwrite = FALSE)
saveRDS(list(implementation = implementation, cohort_dir = normalizePath(cohort_dir),
             annotation_path = normalizePath(annotation_path), module_rule_path = normalizePath(module_rule_path),
             score_members_path = normalizePath(score_members_path), models = models, quantiles = taus,
             h_multiplier = h_multiplier, lambda_gamma = lambda_gamma, random_slope = random_slope,
             subject_ids = unique(as.character(metadata$subject_id)), sample_ids = rownames(metadata)),
        file.path(out_dir, "run_request.rds"))
withCallingHandlers({
all_inference <- list()
for (model in models) {
  n <- length(unique(metadata$subject_id))
  for (tau in taus) {
    inference <- lowdim_profile_pair_workflow_inference_v2(
      metadata, expression_probe, annotation, modules, score_members, module_rule,
      tau = tau, h = h_multiplier * n^(-3 / 10), lambda_gamma = lambda_gamma, include_sledai = model == "B",
      checkpoint_dir = file.path(out_dir, sprintf("workflow_model_%s_tau_%02d", model, round(100 * tau))),
      random_slope = random_slope, checkpoint_storage = "compact"
    )
    write.csv(inference$preprocessing$module_audit,
               file.path(out_dir, sprintf("module_audit_%s_tau_%02d.csv", model, round(100 * tau))),
               row.names = FALSE)
    inference$table$model <- model
    inference$table$tau <- tau
    inference$table$h_multiplier <- h_multiplier
    inference$table <- application_module_simultaneous_intervals_v2(
      inference$table, family_size = if (tau == 0.10) 10L else 30L
    )
    inference$table$simultaneous_family <- if (tau == 0.10) "two_models_tau010_five_modules" else
      "two_models_three_prespecified_quantiles_five_modules"
    all_inference[[length(all_inference) + 1L]] <- inference$table
    artifact_path <- file.path(out_dir, sprintf("model_%s_tau_%02d.rds", model, round(100 * tau)))
    save_compact_application_workflow_v2(inference, artifact_path,
                                         extra = list(implementation = implementation))
    replayed <- read_compact_application_workflow_v2(artifact_path)
    if (!isTRUE(all.equal(replayed$table, inference$table, tolerance = 0, check.attributes = TRUE))) {
      stop("The observed complete workflow table changed after compact readback.")
    }
    combined <- do.call(rbind, all_inference)
    write.csv(combined[combined$coordinate_type == "module", ],
               file.path(out_dir, "confirmatory_module_inference.csv"), row.names = FALSE)
    write.csv(combined[combined$coordinate_type == "clinical", ],
               file.path(out_dir, "clinical_coordinate_inference.csv"), row.names = FALSE)
  }
}
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
jsonlite::write_json(list(
  run_id = run_id, implementation_commit = implementation$commit_sha,
  numeric_certificate_commit = implementation$numeric_certificate_commit,
  actual_source_commit = implementation$actual_source_commit,
  solver_configuration = implementation$solver_configuration,
  cohort_dir = normalizePath(cohort_dir), score_source = score_source,
  method_version = "profile_pair_20261022",
  jackknife_scope = "complete_sample_standardised_workflow",
  parameterisation = "sample_standardised_modules_and_sample_prepared_clinical_design",
  target_type = "unsmoothed_regularised_profile", models = models,
  quantiles = taus, h_multiplier = h_multiplier, lambda_gamma = lambda_gamma,
  random_slope = random_slope,
  simultaneous_module_family_size = sort(unique(combined$simultaneous_family_size[combined$coordinate_type == "module"])),
  completed_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
), file.path(out_dir, "completed_manifest.json"), auto_unbox = TRUE, pretty = TRUE)
message("Low-dimensional profile output: ", normalizePath(out_dir))
}, error = function(condition) {
  saveRDS(condition, file.path(out_dir, "failure_condition.rds"))
  jsonlite::write_json(list(status = "failed", failure_stage = condition$failure_stage %||% "module_workflow",
                            message = conditionMessage(condition), implementation = implementation),
                       file.path(out_dir, "failed_manifest.json"), auto_unbox = TRUE, pretty = TRUE)
})
