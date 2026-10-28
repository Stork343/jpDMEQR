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
input_directory <- file.path(root, "data", "derived", "GSE65391", "subject_design_20261022_v1")
roles_file <- file.path(root, "data", "derived", "GSE65391", "execution_20261022_v1", "exploratory_subject_roles.csv")
inputs_file <- file.path(input_directory, "analysis_inputs.rds")
annotation_file <- file.path(root, "data", "interim", "GSE65391", "probe_annotation_20261022.csv")
inputs <- readRDS(inputs_file)
metadata <- inputs$eligible_metadata
roles <- utils::read.csv(roles_file, stringsAsFactors = FALSE)
audit <- validate_exploratory_subject_roles_v2(metadata, roles)
jobs <- expand.grid(model = c("A", "B"), tau = c(.25, .5, .75),
                     role_id = audit$role_id, stringsAsFactors = FALSE)
jobs <- merge(jobs, audit[c("role_id", "split", "swapped")], by = "role_id", sort = FALSE)
jobs <- jobs[order(jobs$split, jobs$swapped, jobs$model, jobs$tau), ]
rownames(jobs) <- NULL
jobs$job_id <- seq_len(nrow(jobs))
commit <- system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"), stdout = TRUE)
run_id <- get_arg("run-id", "exploratory_pair_20261022")
output <- file.path(root, "results", "raw", "application", run_id)
cache <- file.path(output, "preprocessing")
prepare_only <- identical(get_arg("prepare-only", "false"), "true")
summarise_only <- identical(get_arg("summarise-only", "false"), "true")
checks_path <- get_arg("checks-manifest", file.path(root, "results", "preflight", "application_checks.json"))
implementation <- if (summarise_only || (prepare_only && !file.exists(checks_path))) NULL else
  validate_application_paired_gate_20261022(root, checks_path, inference = TRUE,
                                             method_version = "profile_pair_20261022")
source_inputs <- list(analysis_inputs = normalizePath(inputs_file), expression = inputs$expression_file,
                      annotation = normalizePath(annotation_file), subject_roles = normalizePath(roles_file))
request <- list(numeric_certificate_commit = implementation$numeric_certificate_commit,
                 preparation_source_commit = if (is.null(implementation)) commit else NULL,
                 solver_configuration = implementation$solver_configuration,
                 method_version = "profile_pair_20261022",
                 score_definition = "profile_bandwidth_pair_2_half_minus_full", source_inputs = source_inputs,
                 roles = roles, jobs = jobs, filter_size = 5000L, selected_genes = 3L,
                 mean_score_multiplier = 1, lambda_multiplier = 1, lambda_gamma = 1,
                 ranking_rule = "descending_abs_beta_l1_then_abs_equal_subject_quantile_score_then_gene_id",
                 target_scope = "exploratory_role_specific_selected_design", ci_level = .95)
if (dir.exists(output)) {
  previous_request <- readRDS(file.path(output, "run_request.rds"))
  if (summarise_only) {
    request$numeric_certificate_commit <- previous_request$numeric_certificate_commit
    request$preparation_source_commit <- previous_request$preparation_source_commit
    request$solver_configuration <- previous_request$solver_configuration
  }
  if (!identical(previous_request, request)) stop("The exploratory run request changed.")
} else {
  dir.create(cache, recursive = TRUE)
  saveRDS(request, file.path(output, "run_request.rds"))
  utils::write.csv(jobs, file.path(output, "scheduled_jobs.csv"), row.names = FALSE)
  utils::write.csv(audit, file.path(output, "subject_role_audit.csv"), row.names = FALSE)
}
summarise_run <- function() {
  results <- lapply(jobs$job_id, function(id) {
    path <- file.path(output, sprintf("job%03d", id), "result.rds")
    if (file.exists(path)) readRDS(path) else NULL
  })
  summary <- summarise_exploratory_gene_roles_v2(jobs, results)
  names <- c(role_metrics = "exploratory_role_metrics.csv", intervals = "exploratory_role_intervals.csv",
             selection_frequency = "exploratory_selection_frequency.csv",
             coefficient_median_iqr = "exploratory_coefficient_median_iqr.csv")
  for (name in names(names)) {
    utils::write.csv(summary[[name]], file.path(output, names[[name]]), row.names = FALSE)
  }
  saveRDS(summary$selection, file.path(output, "exploratory_selection_records.rds"))
  invisible(summary)
}
if (summarise_only) {
  summarise_run()
  quit(save = "no", status = 0L)
}
if (prepare_only) {
  expression <- readRDS(inputs$expression_file)[, rownames(metadata), drop = FALSE]
  annotation <- utils::read.csv(annotation_file, stringsAsFactors = FALSE)
  for (i in seq_len(nrow(audit))) {
    path <- file.path(cache, paste0(audit$role_id[i], ".rds"))
    if (file.exists(path)) next
    assignment <- roles[roles$split == audit$split[i] & roles$swapped == audit$swapped[i], ]
    prepared <- prepare_exploratory_gene_split_v2(metadata, expression, annotation, assignment)
    saveRDS(prepared, path)
    message("Prepared ", audit$role_id[i], ": ", ncol(prepared$genes), " genes; clinical A/B ",
            ncol(prepared$clinical$A$X), "/", ncol(prepared$clinical$B$X))
  }
  quit(save = "no", status = 0L)
}
index <- as.integer(get_arg("job-index", "1"))
if (length(index) != 1L || !is.finite(index) || !index %in% jobs$job_id) stop("Supply a declared exploratory job index.")
job <- jobs[jobs$job_id == index, ]
prepared_file <- file.path(cache, paste0(job$role_id, ".rds"))
if (!file.exists(prepared_file)) stop("The complete outcome-blind preprocessing cache is required.")
checkpoint <- file.path(output, sprintf("job%03d", index))
if (dir.exists(checkpoint)) stop("The exploratory job has already been attempted.")
started <- proc.time()[["elapsed"]]
withCallingHandlers({
  result <- fit_exploratory_gene_role_v2(
    metadata, readRDS(prepared_file), job$model, job$tau, checkpoint,
    c(source_inputs, list(preprocessing_file = normalizePath(prepared_file),
                          implementation = implementation))
  )
  jsonlite::write_json(list(job = job, status = result$status,
                            runtime_sec = proc.time()[["elapsed"]] - started,
                            actual_source_commit = commit,
                            implementation = implementation,
                            completed_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)),
                       file.path(checkpoint, "completed_manifest.json"),
                       auto_unbox = TRUE, pretty = TRUE, digits = NA)
}, error = function(condition) {
  jsonlite::write_json(list(job = job, status = "failed", message = conditionMessage(condition),
                            failure_stage = condition$failure_stage %||% "exploratory_fit",
                            runtime_sec = proc.time()[["elapsed"]] - started,
                            actual_source_commit = commit,
                            implementation = implementation), file.path(checkpoint, "failed_manifest.json"),
                       auto_unbox = TRUE, pretty = TRUE, digits = NA)
})
