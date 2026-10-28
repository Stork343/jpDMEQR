root <- normalizePath(".", mustWork = TRUE)
.libPaths(c(file.path(root, "results", "scratch", "profile_pair_package_20261022", "library"), .libPaths()))
library(jpDMEQR)
Sys.setenv(JPDMEQR_CLUSTER_SOLVER = "compiled")
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", envir = environment())
source_v2_module(root, "simulation_v2", envir = environment())
source_v2_module(root, "application_v2", envir = environment())
profile_pair_population_eigen_init_v2(root)
options(jpDMEQR.population_nuisance_backend = "compiled")
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}
implementation <- validate_application_population_implementation_v2(root, get_arg("core-build-commit", "5c47a62"))
config <- yaml::read_yaml(file.path(root, "config", "application_matched_prospective_20261022.yml"))
model <- get_arg("model", "A")
q <- as.integer(get_arg("q", "1"))
error <- get_arg("error", "normal")
if (!model %in% config$models || !q %in% config$random_effect_dimensions ||
    !error %in% config$error_families) stop("Unknown matched population case.")
case <- sprintf("AM_%s_q%d_%s", model, q, error)
case_index <- match(model, config$models) + 2L * (q - 1L) + 4L * (match(error, config$error_families) - 1L)
population <- readRDS(file.path(root, "data", "derived", "GSE65391", "execution_20261022_v1",
                                 sprintf("matched_population_model_%s_q%d.rds", model, q)))
run_id <- get_arg("run-id", paste0(case, "_targets_", format(Sys.time(), "%Y%m%dT%H%M%S")))
output <- file.path(root, "results", "raw", "application_matched", run_id)
if (dir.exists(output)) stop("The matched target output already exists.")
dir.create(output, recursive = TRUE)
n_population <- config$population_targets$clusters_per_repeat
repeats <- config$population_targets$independent_repeats
h_analysis <- config$n_subjects^(-3 / 10)
h_target <- min(h_analysis / 4, n_population^(-1 / 3))
control <- profile_pair_asset_control_v2(list(population_engine = get_arg("engine", "streamed")))
estimates <- refinements <- matrix(NA_real_, repeats, length(population$beta0),
                                    dimnames = list(NULL, names(population$beta0)))
records <- vector("list", repeats)
for (repeat_id in seq_len(repeats)) {
  seed <- config$base_seed + 100000000L + case_index * 1000000L + repeat_id * 10000L
  dat <- generate_application_matched_population_data_v2(population, config, error, n_population, seed)
  if (control$population_engine == "streamed") dat$population_layout <- profile_pair_population_layout_v2(dat)
  lower <- solve_profile_pair_target_v2(dat, 0.25, h_target, 1,
                                       beta_start = population$beta0, control = control)
  saveRDS(profile_pair_fit_record_v2(lower), file.path(output, sprintf("lower_repeat_%02d.rds", repeat_id)))
  refinement <- solve_profile_pair_target_v2(dat, 0.25, h_target / 2, 1,
                                            beta_start = lower$beta, control = control)
  saveRDS(profile_pair_fit_record_v2(refinement), file.path(output, sprintf("refinement_repeat_%02d.rds", repeat_id)))
  estimates[repeat_id, ] <- lower$beta
  refinements[repeat_id, ] <- refinement$beta
  records[[repeat_id]] <- list(seed = seed, n_population = n_population,
                               lower = profile_pair_fit_record_v2(lower),
                               refinement = profile_pair_fit_record_v2(refinement))
  rm(dat, lower, refinement)
  invisible(gc())
}
lower_target <- colMeans(estimates)
upper_target <- 2 * population$beta0 - lower_target
audits <- list()
for (repeat_id in seq_len(repeats)) {
  seed <- config$base_seed + 200000000L + case_index * 1000000L + repeat_id * 10000L
  dat <- generate_application_matched_population_data_v2(population, config, error, n_population, seed)
  if (control$population_engine == "streamed") dat$population_layout <- profile_pair_population_layout_v2(dat)
  for (tau in c(0.5, 0.75)) {
    beta <- if (tau == 0.5) population$beta0 else upper_target
    components <- profile_pair_population_evaluate_v2(dat, beta, tau, h_target, 1, control, 2L)
    moments <- profile_pair_moments_v2(components)
    audits[[length(audits) + 1L]] <- data.frame(
      repeat_id = repeat_id, seed = seed, tau = tau, coordinate = names(beta),
      score_mean = as.numeric(moments$score_mean),
      score_mcse = sqrt(diag(moments$Sigma_population) / n_population),
      nuisance_gradient = components$max_nuisance_gradient,
      profile_identity_error = components$profile_identity_error,
      stringsAsFactors = FALSE
    )
  }
  rm(dat, components)
  invisible(gc())
}
audit <- do.call(rbind, audits)
mcse <- apply(estimates, 2, stats::sd) / sqrt(repeats)
target <- rbind(
  data.frame(tau = .25, coordinate = names(lower_target), target = as.numeric(lower_target),
              target_mcse = mcse, target_source = "independent_complete_profile_mc"),
  data.frame(tau = .50, coordinate = names(population$beta0), target = as.numeric(population$beta0),
              target_mcse = 0, target_source = "analytic_central_symmetry"),
  data.frame(tau = .75, coordinate = names(upper_target), target = as.numeric(upper_target),
              target_mcse = mcse, target_source = "derived_exact_reflection_of_lower_target")
)
target$case_id <- case
write.csv(target, file.path(output, "population_targets.csv"), row.names = FALSE)
write.csv(audit, file.path(output, "independent_population_score_audit.csv"), row.names = FALSE)
saveRDS(list(case = case, population_beta = population$beta0, targets = target, estimates = estimates,
             refinements = refinements, target_records = records, score_audit = audit,
             population_reference = population[c("X", "Z", "rows", "beta0", "preprocessing")],
             method_version = "profile_pair_20261022", implementation_commit = implementation$commit_sha,
             core_build_commit = implementation$core_build_commit, population_core_checks = implementation,
             h_target = h_target, h_analysis = h_analysis, config = config),
        file.path(output, "matched_target_asset.rds"))
jsonlite::write_json(list(
  case = case, method_version = "profile_pair_20261022", implementation_commit = implementation$commit_sha,
  core_build_commit = implementation$core_build_commit,
  n_population = n_population, independent_repeats = repeats, h_target = h_target,
  h_refinement = h_target / 2, target_mcse_max = max(mcse),
  target_refinement_max = max(abs(refinements - estimates)),
  population_engine = control$population_engine,
  completed_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
), file.path(output, "completed_manifest.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
