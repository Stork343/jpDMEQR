root <- normalizePath(".", mustWork = TRUE)
.libPaths(c(normalizePath(file.path(root, "results", "scratch", "profile_pair_package_20261022", "library")), .libPaths()))
library(jpDMEQR)
source(file.path(root, "scripts", "00_source_v2.R"))
for (module in c("profile_v2", "simulation_v2", "metrics_v2")) source_v2_module(root, module)
source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"))
cli <- parse_cli_args_v2(commandArgs(trailingOnly = TRUE))
commit <- profile_pair_asset_commit_v2(cli$commit %||% stop("Supply --commit."), FALSE)
if (!profile_pair_population_commit_matches_v2(commit, current_commit_v2(root), "streamed", root)) {
  stop("Structural-audit population dependencies differ from the declared commit.")
}
run_id <- cli$`run-id` %||% stop("Supply a new structural-audit run-id.")
workers <- as.integer(cli$workers %||% "2")
if (!is.finite(workers) || workers < 1L || workers > 2L) stop("Use one or two audit workers.")
out <- file.path(root, "results", "scratch", "profile_pair_structural_audits_20261022", run_id)
if (dir.exists(out)) stop("The structural-audit run directory already exists.")
dir.create(out, recursive = TRUE)
registry <- read_simulation_registry_v2(file.path(root, "config", "simulation_pair_main_20261022.yml"))
configs <- lapply(c("E03", "E04", "E05", "E06"), function(id) {
  config_row_to_list_v2(registry[registry$experiment_id == id, , drop = FALSE])
})
if (!setequal(vapply(configs, `[[`, numeric(1), "lambda_gamma"), c(0.25, 0.5, 2, 4))) {
  stop("The frozen nuisance-penalty audit grid is incomplete.")
}
baseline <- configs[[1]]
baseline$experiment_id <- "LAMBDA1_BASELINE"
baseline$lambda_gamma <- 1
configs <- c(configs, list(baseline))
Sys.setenv(JPDMEQR_CLUSTER_SOLVER = "newton", JPDMEQR_DANTZIG_SOLVER = "clarabel_native",
           JPDMEQR_ASSEMBLY = "sequential", JPDMEQR_CV_CORES = "1",
           JPDMEQR_HESSIAN_REFRESH = "1", JPDMEQR_SHARED_FIRST_STAGE = "FALSE",
           OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
profile_pair_population_eigen_init_v2(root)
options(jpDMEQR.population_nuisance_backend = "compiled")
audit_dir <- file.path(root, "results", "intermediate", "profile_pair_20261022", "construction_v3", "structural_audits")
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(configs, file.path(out, "audit_designs.rds"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
audit_one <- function(registered_config) {
  id <- registered_config$experiment_id
  path <- file.path(audit_dir, paste0(id, ".rds"))
  if (file.exists(path)) stop("A structural audit already exists: ", path)
  case_dir <- file.path(out, id)
  dir.create(case_dir, recursive = TRUE)
  stage <- "low_dimensional_audit_design"
  withCallingHandlers({
    source_config <- if (identical(id, "LAMBDA1_BASELINE")) {
      config_row_to_list_v2(registry[registry$experiment_id == "E03", , drop = FALSE])
    } else registered_config
    columns <- population_target_columns_v2(registered_config)
    if (!identical(columns, seq_len(max(columns)))) stop("Audit columns must form the declared initial marginal design.")
    config <- registered_config
    config$p <- max(columns)
    stage <- "paired_population_roots"
    checkpoint <- function(label, object) {
      stage <<- label
      saveRDS(object, file.path(case_dir, paste0(label, ".rds")), compress = FALSE)
    }
    target <- build_profile_pair_target_asset_v2(
      config, n_population = 100000L, repeats = 4L, development = FALSE,
      implementation_commit = commit,
      control = list(population_engine = "streamed", batch_clusters = 1000L), checkpoint = checkpoint)
    beta0 <- make_beta_v2(config$p, config$s, config$signal)
    object <- list(
      asset_kind = "paired_structural_target_audit", audit_construction_version = "paired-structural-audit-20261022-v1",
      method_version = target$method_version, score_definition = target$score_definition,
      construction_version = target$construction_version,
      implementation_commit = commit, development = FALSE,
      experiment_id = id, registered_p = registered_config$p, audit_p = config$p,
      audit_columns = columns, registered_config = source_config, audit_design = registered_config,
      registered_experiment_id = source_config$experiment_id,
      audit_design_origin = if (identical(id, "LAMBDA1_BASELINE")) {
        "frozen_E_design_with_requested_Lambda1_reference"
      } else "registered_E_sensitivity",
      config = config,
      audit_scope = "all_active_coefficients_and_three_null_controls",
      beta0 = beta0, beta_star_mc = target$beta_star_mc,
      target_displacement = target$beta_star_mc - beta0,
      target_mc_se_by_coordinate = target$target_mc_se_by_coordinate,
      target_mc_se = target$target_mc_se, estimates_by_repeat = target$estimates_by_repeat,
      target_displacement_by_repeat = sweep(target$estimates_by_repeat, 2, beta0),
      between_run_max_difference = target$between_run_max_difference,
      n_population = target$n_population, repeats = target$repeats,
      all_converged = target$all_converged, max_kkt_residual = target$max_kkt_residual,
      max_nuisance_gradient = target$max_nuisance_gradient,
      profile_identity_error = target$profile_identity_error,
      h_target = target$h_target, h_analysis = target$h_analysis,
      paired_target_asset = target,
      created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
    validate_profile_pair_target_asset_v2(target, config, final = TRUE, expected_commit = current_commit_v2(root))
    check <- validate_structural_target_audit_v2(object, final = TRUE)
    saveRDS(object, path, compress = FALSE)
    if (!check$valid) stop(paste(check$problems, collapse = "; "))
    restored <- readRDS(path)
    validate_profile_pair_target_asset_v2(restored$paired_target_asset, config, final = TRUE,
                                          expected_commit = current_commit_v2(root))
    if (!validate_structural_target_audit_v2(restored, final = TRUE)$valid) stop("Structural-audit read-back failed.")
    coordinate <- data.frame(
      experiment_id = id, lambda_gamma = registered_config$lambda_gamma,
      coordinate = names(beta0), beta0 = beta0, beta_star_mc = target$beta_star_mc,
      displacement = object$target_displacement, mcse = target$target_mc_se_by_coordinate)
    utils::write.csv(coordinate, file.path(case_dir, "coordinates.csv"), row.names = FALSE)
    jsonlite::write_json(list(pass = TRUE, path = path, experiment_id = id,
                              lambda_gamma = registered_config$lambda_gamma,
                              registered_p = registered_config$p, audit_p = config$p,
                              n_population = 100000L, repeats = 4L,
                              max_absolute_displacement = max(abs(object$target_displacement)),
                              max_target_mc_se = target$target_mc_se,
                              kkt = target$max_kkt_residual,
                              nuisance = target$max_nuisance_gradient),
                         file.path(case_dir, "manifest.json"), pretty = TRUE, auto_unbox = TRUE, digits = NA)
    message("Structural audit completed: ", id)
    list(experiment_id = id, path = path, valid = TRUE)
  }, error = function(error) {
    saveRDS(list(condition = error, stage = stage), file.path(case_dir, "failure.rds"))
    writeLines(c(stage, conditionMessage(error)), file.path(case_dir, "FAILED.txt"))
  })
}
result <- if (workers == 1L) lapply(configs, audit_one) else {
  parallel::mclapply(configs, audit_one, mc.cores = workers, mc.preschedule = FALSE)
}
saveRDS(result, file.path(out, "build_results.rds"))
if (any(vapply(result, inherits, logical(1), "try-error"))) stop("A structural audit failed; inspect its saved condition.")
jsonlite::write_json(list(pass = TRUE, implementation_commit = commit, cases = result,
                          n_population = 100000L, repeats = 4L, workers = workers),
                     file.path(out, "manifest.json"), pretty = TRUE, auto_unbox = TRUE, digits = NA)
writeLines("COMPLETE", file.path(out, "COMPLETE.txt"))
