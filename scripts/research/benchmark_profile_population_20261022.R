root <- normalizePath(".", mustWork = TRUE)
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2")
source_v2_module(root, "simulation_v2")
Sys.setenv(JPDMEQR_CLUSTER_SOLVER = "newton")
backend <- Sys.getenv("JPDMEQR_POPULATION_TEST_BACKEND", unset = "R_BLAS")
if (backend %in% c("Eigen", "compiled_Eigen")) {
  .libPaths(c(normalizePath(file.path(root, "results", "scratch", "profile_pair_package_20261022", "library")), .libPaths()))
  library(jpDMEQR)
  profile_pair_population_eigen_init_v2(root)
  if (backend == "compiled_Eigen") options(jpDMEQR.population_nuisance_backend = "compiled")
}
out <- file.path(root, "results", "scratch", "profile_pair_population_performance_20261022")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- file.path(out, backend)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
registry <- read_simulation_registry_v2(file.path(root, "config", "simulation_main.csv"))
rows <- list()
for (id in c("C14", "A05", "A08")) {
  config <- config_row_to_list_v2(registry[registry$experiment_id == id, , drop = FALSE])
  n <- if (config$p == 300L) 2000L else if (config$p == 1000L) 1000L else 500L
  message("Actual full-design component benchmark: ", id, " p=", config$p, " n=", n)
  started <- proc.time()[[3]]
  dat <- profile_pair_population_data_v2(config, n, 20261211L)
  generation <- proc.time()[[3]] - started
  dat$population_layout <- profile_pair_population_layout_v2(dat)
  started <- proc.time()[[3]]
  result <- profile_pair_population_streamed_v2(
    dat, dat$beta0, config$tau, config$n_clusters^(-3 / 10), config$lambda_gamma,
    batch_clusters = 250L)
  elapsed <- proc.time()[[3]] - started
  rows[[length(rows) + 1L]] <- data.frame(
    experiment_id = id, p = config$p, n_clusters = n,
    generation_sec = generation, component_sec = elapsed,
    linear_scale_100000_sec = elapsed * 100000 / n,
    observed_data_bytes = as.numeric(object.size(dat)),
    full_100000_X_bytes = 100000 * mean(config$m_values) * config$p * 8,
    max_nuisance_gradient = result$max_nuisance_gradient,
    identity_error = result$profile_identity_error,
    covariance_identity_error = result$population_moments$covariance_identity_error)
  utils::write.csv(do.call(rbind, rows), file.path(out, "full_dimension_timing.csv"), row.names = FALSE)
  rm(dat, result)
  gc()
}
comparisons <- list()
base <- config_row_to_list_v2(registry[registry$experiment_id == "D01", , drop = FALSE])
base$p <- 12L
base$s <- 3L
for (design in c("ar1", "independent", "block", "sparse_precision", "factor", "dense_precision")) {
  config <- base
  config$design_type <- design
  for (rr in 1:4) {
    seed <- 20261220L + rr
    dat <- profile_pair_population_data_v2(config, 2000L, seed)
    integrated <- profile_pair_conditional_gaussian_data_v2(config, 2000L, seed)
    finite <- profile_pair_population_streamed_v2(dat, dat$beta0, 0.5, 0.2, 1,
                                                  batch_clusters = 500L)
    conditional <- profile_pair_conditional_gaussian_components_v2(integrated, integrated$beta0, 0.5, 0.2, 1)
    comparisons[[length(comparisons) + 1L]] <- data.frame(
      design = design, repeat_id = rr, n_clusters = 2000L, p = 12L,
      nuisance_observation_error = max(abs(integrated$residual_base - (dat$y - dat$X %*% dat$beta0))),
      H11_difference = finite$hessian[1, 1] - conditional$hessian[1, 1],
      H12_difference = finite$hessian[1, 2] - conditional$hessian[1, 2],
      covariance11_difference = finite$population_moments$Sigma_population[1, 1] -
        conditional$population_moments$Sigma_population[1, 1],
      finite_score1 = finite$score[1], conditional_score1 = conditional$score[1])
  }
  message("Conditional-versus-finite comparison completed: ", design)
}
comparison <- do.call(rbind, comparisons)
utils::write.csv(comparison, file.path(out, "conditional_finite_MC_comparison.csv"), row.names = FALSE)
summaries <- lapply(split(comparison, comparison$design), function(d) {
  fields <- c("H11_difference", "H12_difference", "covariance11_difference", "finite_score1")
  do.call(rbind, lapply(fields, function(field) data.frame(
    design = d$design[1], quantity = field, repeats = nrow(d),
    average = mean(d[[field]]), mcse = stats::sd(d[[field]]) / sqrt(nrow(d)))))
})
utils::write.csv(do.call(rbind, summaries), file.path(out, "conditional_finite_MC_summary.csv"), row.names = FALSE)
writeLines("COMPLETE", file.path(out, "benchmark_COMPLETE.txt"))
print(do.call(rbind, rows))
