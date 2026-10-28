root <- normalizePath(".", mustWork = TRUE)
.libPaths(c(file.path(root, "results", "scratch", "profile_pair_package_20261022", "library"), .libPaths()))
library(jpDMEQR)
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", environment())
source_v2_module(root, "application_v2", environment())
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Supply the declared case ID and a new equivalence run ID.")
registry <- read.csv(file.path(root, "data", "derived", "GSE65391", "execution_20261022_v1", "matched_registry.csv"),
                      stringsAsFactors = FALSE)
case <- registry[registry$case_id == args[1], ]
if (nrow(case) != 1L || !case$case_id %in% c("AM_A_q1_normal_tau25", "AM_B_q2_t3_tau75")) {
  stop("Use one of the two declared complete-workflow equivalence cases.")
}
output <- file.path(root, "results", "scratch", args[2])
if (dir.exists(output)) stop("The complete configuration equivalence directory already exists.")
dir.create(output, recursive = TRUE)
original_path <- file.path(root, "results", "scratch", "numerical_solver_acceptance_20261022_v3", "solver_configuration.json")
production_path <- file.path(root, "results", "scratch", "profile_pair_geometry_development_20261022",
                              "production_validation_20261022_02", "manifest.json")
configurations <- list(
  original148 = jsonlite::read_json(original_path, simplifyVector = TRUE),
  production = jsonlite::read_json(production_path, simplifyVector = TRUE)$solver_configuration
)
set_configuration <- function(configuration) {
  environment_names <- c(
    dantzig_solver = "JPDMEQR_DANTZIG_SOLVER", cluster_solver = "JPDMEQR_CLUSTER_SOLVER",
    assembly = "JPDMEQR_ASSEMBLY", cv_workers = "JPDMEQR_CV_CORES",
    hessian_refresh = "JPDMEQR_HESSIAN_REFRESH", unpenalised_optimizer = "JPDMEQR_UNPENALISED_OPTIMIZER",
    proximal_scaling = "JPDMEQR_PROXIMAL_SCALING", shared_first_stage = "JPDMEQR_SHARED_FIRST_STAGE",
    matrix_backend = "JPDMEQR_MATRIX_BACKEND", identity_detail = "JPDMEQR_IDENTITY_DETAIL"
  )
  values <- lapply(configuration[names(environment_names)], as.character)
  names(values) <- unname(environment_names)
  do.call(Sys.setenv, values)
  options(jpDMEQR.matrix_backend = configuration$matrix_backend,
          jpDMEQR.identity_detail = configuration$identity_detail)
  if (!configuration_matches_v2(configuration)) stop("The actual requested configuration was not installed.")
}
set_configuration(configurations$original148)
source_identity <- validate_application_numeric_source_identity_v2(root, "078559a")
config <- yaml::read_yaml(file.path(root, "config", "application_matched_prospective_20261022.yml"))
population <- readRDS(file.path(root, "data", "derived", "GSE65391", "execution_20261022_v1",
                                 sprintf("matched_population_model_%s_q%d.rds", case$model, case$q)))
rule <- yaml::read_yaml(file.path(root, "config", "module_scoring_20261022.yml"))
members <- read_frozen_module_score_members_v2(file.path(root, rule$score_membership_file))
signatures <- read_frozen_module_definitions_v2(file.path(root, rule$membership_file))
annotation <- read.csv(file.path(root, rule$annotation_path), stringsAsFactors = FALSE)
expression <- readRDS(file.path(root, "data", "derived", "GSE65391", "module_measurements_20261022_v1",
                                 "observed_module_probe_expression.rds"))
seed <- config$base_seed + (case$q + 2L * (match(case$error_family, config$error_families) - 1L)) * 1000000L + 1L
draw <- resample_application_matched_workflow_v2(population, expression, config, case$error_family, seed)
saveRDS(draw, file.path(output, "common_observed_expression_and_declared_draw.rds"))
saveRDS(list(case = case, seed = seed, source_identity = source_identity,
             configuration_sources = c(original148 = original_path, production = production_path)),
        file.path(output, "validation_request.rds"))
asset <- readRDS(file.path(root, "results", "raw", "application_matched",
                           sprintf("matched_targets_20261022_v1_AM_%s_q%d_%s", case$model, case$q, case$error_family),
                           "matched_target_asset.rds"))
target <- asset$targets[asset$targets$tau == case$tau, ]
results <- mapped <- list()
configuration_summary <- list()
for (name in names(configurations)) {
  set_configuration(configurations[[name]])
  directory <- file.path(output, name)
  dir.create(directory)
  jsonlite::write_json(solver_configuration_v2(), file.path(directory, "actual_configuration.json"),
                       auto_unbox = TRUE, pretty = TRUE)
  started <- proc.time()[[3]]
  result <- lowdim_profile_pair_workflow_inference_v2(
    draw$metadata, draw$expression, annotation, signatures, members, rule,
    case$tau, 129^(-3 / 10), include_sledai = case$model == "B", random_slope = case$q == 2L,
    checkpoint_dir = file.path(directory, "workflow_checkpoints"), checkpoint_storage = "compact"
  )
  elapsed <- proc.time()[[3]] - started
  mapped[[name]] <- application_matched_coordinate_metrics_v2(result, population, target, case, 1L)
  save_compact_application_workflow_v2(result, file.path(directory, "workflow_result.rds"),
                                       extra = list(mapped_jackknife = mapped[[name]]$mapped_jackknife,
                                                    actual_configuration = solver_configuration_v2()))
  results[[name]] <- read_compact_application_workflow_v2(file.path(directory, "workflow_result.rds"))
  diagnostics <- rbind(result$jackknife$full_sample$diagnostics, result$jackknife$diagnostics)
  write.csv(diagnostics, file.path(directory, "all_fit_diagnostics.csv"), row.names = FALSE)
  configuration_summary[[name]] <- data.frame(
    configuration = name, n_subjects = 129L, n_visits = nrow(draw$metadata), p = ncol(population$X),
    completed_fits = nrow(diagnostics), completed_jackknives = result$jackknife$n_clusters,
    runtime_sec = elapsed, max_refit_gradient = max(diagnostics$refit_gradient),
    max_nuisance_gradient = max(diagnostics$nuisance_gradient),
    max_identity_error = max(diagnostics$profile_identity_error), seed = seed,
    actual_source_commit = source_identity$actual_source_commit
  )
}
comparisons <- list()
compare <- function(name, reference, candidate, tolerance, relative = FALSE) {
  if (!identical(dim(reference), dim(candidate)) || length(reference) != length(candidate) ||
      any(!is.finite(reference)) || any(!is.finite(candidate))) stop("Invalid equivalence comparison: ", name)
  difference <- max(abs(reference - candidate))
  scale <- if (relative) max(1, max(abs(reference))) else 1
  comparisons[[length(comparisons) + 1L]] <<- data.frame(
    quantity = name, max_absolute_difference = difference, comparison_scale = scale,
    tolerance = tolerance, passed = difference <= tolerance * scale
  )
}
a <- results$original148
b <- results$production
compare("beta_refit", a$fit$beta_refit, b$fit$beta_refit, 1e-8, TRUE)
compare("beta_pair", a$fit$beta_pair, b$fit$beta_pair, 1e-8, TRUE)
compare("fitted_paired_score", a$fit$components$score, b$fit$components$score, 2e-4)
compare("fitted_paired_Hessian", a$fit$components$hessian, b$fit$components$hessian, 3e-3)
compare("raw_jackknife_matrix", a$jackknife$estimates, b$jackknife$estimates, 1e-8, TRUE)
compare("population_unit_jackknife_matrix", mapped$original148$mapped_jackknife,
        mapped$production$mapped_jackknife, 1e-8, TRUE)
pointwise <- c("beta_hat", "sandwich_se", "sandwich_lower", "sandwich_upper", "jackknife_se", "jackknife_lower", "jackknife_upper")
compare("all_pointwise_intervals", as.matrix(mapped$original148$table[pointwise]),
        as.matrix(mapped$production$table[pointwise]), 1e-8, TRUE)
family <- c("sandwich_simultaneous_lower", "sandwich_simultaneous_upper", "jackknife_simultaneous_lower", "jackknife_simultaneous_upper")
module <- grepl("^module__", mapped$original148$table$coordinate)
compare("family30_intervals", as.matrix(mapped$original148$table[module, family]),
        as.matrix(mapped$production$table[module, family]), 1e-8, TRUE)
same_preprocessing <- identical(a$preprocessing, b$preprocessing) &&
  all(vapply(seq_along(a$jackknife$leave_one), function(i) {
    identical(a$jackknife$leave_one[[i]]$preprocessing, b$jackknife$leave_one[[i]]$preprocessing)
  }, logical(1)))
if (!same_preprocessing) stop("A configuration changed one of the complete 130 preprocessing definitions.")
prepared <- prepare_gse65391_profile_workflow_stage_v2(
  draw$metadata, draw$expression, annotation, signatures, members, rule,
  include_sledai = case$model == "B", random_slope = case$q == 2L
)
common <- list()
for (name in names(configurations)) {
  set_configuration(configurations[[name]])
  common[[name]] <- profile_bandwidth_pair_v2(
    prepared$y, prepared$X, prepared$Z, prepared$cluster_id, a$fit$beta_refit, case$tau,
    129^(-3 / 10), lambda_gamma = 1, nuisance_control = list(grad_tol = 1e-11), identity_check = TRUE
  )
}
compare("common_beta_paired_score", common$original148$score, common$production$score, 1e-7)
compare("common_beta_cluster_scores", common$original148$cluster_scores, common$production$cluster_scores, 1e-7)
compare("common_beta_Hessian", common$original148$hessian, common$production$hessian, 1e-8)
table <- do.call(rbind, comparisons)
summary <- do.call(rbind, configuration_summary)
write.csv(table, file.path(output, "equivalence_comparisons.csv"), row.names = FALSE)
write.csv(summary, file.path(output, "configuration_fit_summary.csv"), row.names = FALSE)
if (any(!table$passed) || any(summary$completed_fits != 130L) ||
    any(summary$max_refit_gradient > 1e-7 | summary$max_nuisance_gradient > 1e-7 | summary$max_identity_error > 1e-8)) {
  stop("The complete dual-configuration workflow comparison did not meet the existing numerical standards.")
}
jsonlite::write_json(list(case = case, seed = seed, source_identity = source_identity,
                          configurations = configurations, same_complete_preprocessing = same_preprocessing,
                          passed_comparisons = sum(table$passed),
                          completed_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)),
                     file.path(output, "equivalence_manifest.json"), auto_unbox = TRUE, pretty = TRUE)
