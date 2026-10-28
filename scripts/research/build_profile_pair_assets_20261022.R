args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1]] else ".", mustWork = TRUE)
run_id <- if (length(args) >= 2L) args[[2]] else format(Sys.time(), "%Y%m%dT%H%M%S", tz = "UTC")
target_mode <- if (length(args) >= 3L) args[[3]] else "profile_mc"
if (!target_mode %in% c("structural", "profile_mc")) stop("Unknown target_mode.")
run_dir <- file.path(root, "results", "scratch", "profile_pair_asset_development_20261022", run_id)
if (dir.exists(run_dir)) stop("A development run directory already exists: ", run_dir)
dir.create(run_dir, recursive = TRUE)
stage <- "source"
options(error = function() {
  writeLines(c(paste("stage:", stage), geterrmessage()), file.path(run_dir, "FAILED.txt"))
  writeLines(capture.output(traceback(15)), file.path(run_dir, "traceback.txt"))
  quit(save = "no", status = 1L)
})
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2")
source_v2_module(root, "simulation_v2")
Sys.setenv(JPDMEQR_CLUSTER_SOLVER = "newton", JPDMEQR_DANTZIG_SOLVER = "clarabel_native",
           JPDMEQR_ASSEMBLY = "sequential", JPDMEQR_CV_CORES = "1",
           JPDMEQR_HESSIAN_REFRESH = "1", JPDMEQR_SHARED_FIRST_STAGE = "FALSE")
config <- list(
  experiment_id = paste0("PAIR-ASSET-DEV-", run_id), module = "paired_asset_development",
  n_clusters = 80L, p = 2L, s = 1L, tau = 0.25, q = 1L,
  m_values = 3:5, m_rule = "uniform", m_multiset = NULL,
  design_type = "ar1", rho_x = 0.3, signal = 0.75, error_dist = "gaussian",
  heteroskedastic = FALSE, error_dependence = "independent", copula_rho = 0,
  random_effect_dist = "normal", sigma_b0 = 0.8, sigma_b1 = 0.5,
  rho_b = 0.4, x_b_corr = 0, informative_size = FALSE, nonlinear_re_strength = 0,
  fit_random_effects = "intercept", lambda_gamma = 1, h_multiplier = 1,
  target_mode = target_mode
)
commit <- system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"), stdout = TRUE)
if (length(commit) != 1L || !nzchar(commit)) stop("An implementation commit is required.")
checkpoint <- function(label, object) {
  stage <<- label
  saveRDS(object, file.path(run_dir, paste0(label, ".rds")))
}
saveRDS(config, file.path(run_dir, "config_snapshot.rds"))
writeLines(capture.output(sessionInfo()), file.path(run_dir, "sessionInfo.txt"))
stage <- "target_construction"
target <- build_profile_pair_target_asset_v2(
  config, n_population = 120L, repeats = 2L, seed = 20261022L,
  development = TRUE, implementation_commit = commit, checkpoint = checkpoint
)
target_path <- profile_pair_target_asset_path_v2(root, config$experiment_id)
direction_path <- profile_pair_direction_asset_path_v2(root, config$experiment_id)
if (file.exists(target_path) || file.exists(direction_path)) stop("A paired asset path already exists.")
dir.create(dirname(target_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(direction_path), recursive = TRUE, showWarnings = FALSE)
saveRDS(target, target_path)
stage <- "direction_construction"
direction <- build_profile_pair_direction_asset_v2(
  config, target, n_population = 120L, repeats = 2L, seed = 20271022L,
  development = TRUE, implementation_commit = commit, checkpoint = checkpoint
)
saveRDS(direction, direction_path)
stage <- "read_back_validation"
target_read <- readRDS(target_path)
direction_read <- readRDS(direction_path)
validate_profile_pair_target_asset_v2(target_read, config, final = FALSE, expected_commit = commit)
validate_profile_pair_direction_asset_v2(direction_read, config, final = FALSE, expected_commit = commit)
diagnostic_data <- profile_pair_population_data_v2(config, 40L, seed = 20281022L)
diagnostic <- profile_pair_target_scores_v2(diagnostic_data, target_read)
saveRDS(diagnostic, file.path(run_dir, "independent_target_score_diagnostic.rds"))
diagnostic_direction <- profile_bandwidth_pair_v2(
  diagnostic_data$y, diagnostic_data$X, diagnostic_data$Z, diagnostic_data$cluster_id,
  direction_read$evaluation_point, config$tau, direction_read$h_analysis,
  config$lambda_gamma, nuisance_control = list(grad_tol = 1e-7),
  identity_check = TRUE, need_fold_sums = FALSE
)
saveRDS(diagnostic_direction, file.path(run_dir, "independent_direction_score_diagnostic.rds"))
utils::write.csv(data.frame(
  coordinate = names(target$beta_star_mc), beta_star_mc = target$beta_star_mc,
  target_mc_se = target$target_mc_se_by_coordinate,
  beta_refinement_mc = target$beta_refinement_mc,
  beta_analysis_bandwidth_mc = target$beta_analysis_bandwidth_mc,
  sigma_population = direction$sigma_population
), file.path(run_dir, "coordinate_summary.csv"), row.names = FALSE)
manifest <- list(
  run_id = run_id, experiment_id = config$experiment_id, development = TRUE,
  final_eligible = FALSE, target_path = target_path, direction_path = direction_path,
  score_definition = target$score_definition,
  construction_version = target$construction_version,
  implementation_commit = commit, n_population = target$n_population,
  repeats = target$repeats, h_target = target$h_target, h_analysis = target$h_analysis,
  target_mode = target_mode, direction_target_source = direction$target_source,
  direction_evaluation_point = direction$evaluation_point,
  paired_mc_target = target$beta_star_mc,
  target_max_kkt = target$max_kkt_residual,
  max_nuisance_gradient = max(target$max_nuisance_gradient, direction$max_nuisance_gradient),
  identity_error = max(target$profile_identity_error, direction$profile_identity_error),
  covariance_identity_error = direction$covariance_identity_error,
  target_mc_se = target$target_mc_se,
  target_refinement_distance = target$target_refinement_distance_by_repeat,
  smoothing_target_distance = target$smoothing_target_distance,
  direction_min_eigenvalue = direction$curvature$min_eigenvalue,
  direction_condition = direction$condition_number,
  maximum_inverse_defect = max(vapply(direction$directions, `[[`, numeric(1), "residual")),
  solver_configuration = direction$solver_configuration
)
jsonlite::write_json(manifest, file.path(run_dir, "completed_manifest.json"),
                     pretty = TRUE, auto_unbox = TRUE, digits = NA)
writeLines("COMPLETE", file.path(run_dir, "COMPLETE.txt"))
print(manifest)
