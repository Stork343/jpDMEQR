source("research/submission_completion_20261027/load.R")
reference <- readRDS("results/scratch/submission_comparator_probes_20261027/SQR-DEBIASED-IID_A01_matrix_free_fista/result.rds")$fit
registry <- read_simulation_registry_v2("config/simulation_main.csv")
config <- config_row_to_list_v2(registry[registry$experiment_id == "A01", ])
config$method_version <- "profile_pair_20261022"
seed <- seed_from_id_v2(20261027L, "A01", 1L, "numerical_probe")
data <- prepare_analysis_data_v2(simulate_from_config_v2(config, seed), config, seed)$data
role <- reference$tuning_diagnostics$gamma_cv$roles[[1L]]
training <- subset_clusters_v2(data, role$training_subjects)
initial <- reference$tuning_diagnostics$gamma_cv$fold_initial_estimators[[1L]]
H <- sqr_iid_hessian_fn_v2(training$y, training$X, initial$beta, config$tau, initial$h)
p <- nrow(H)
canonical <- c(1, rep(0, p - 1L))
matrix <- rbind(cbind(H, -H), cbind(-H, H), -diag(2L * p))
output <- "results/scratch/submission_completion_20261027/precision_row_checks"
if (dir.exists(output)) stop("The precision row check already exists.")
dir.create(output, recursive = TRUE)
saveRDS(H, file.path(output, "H.rds"))
records <- list()
for (gamma in unique(reference$tuning_diagnostics$gamma_cv$cells$gamma)[1:6]) {
  started <- proc.time()[[3L]]
  fit <- clarabel::clarabel(A = Matrix::Matrix(matrix, sparse = TRUE),
    b = c(canonical + gamma, -canonical + gamma, rep(0, 2L * p)),
    q = rep(1, 2L * p), cones = list(l = 4L * p),
    control = list(verbose = FALSE, max_iter = 200L, tol_gap_abs = 1e-8,
      tol_gap_rel = 1e-8, tol_feas = 1e-8))
  elapsed <- proc.time()[[3L]] - started
  omega <- fit$x[seq_len(p)] - fit$x[p + seq_len(p)]
  record <- data.frame(gamma = gamma, elapsed = elapsed, status = fit$status,
    residual = max(abs(H %*% omega - canonical)), objective = sum(abs(omega)))
  records[[length(records) + 1L]] <- record
  saveRDS(fit, file.path(output, paste0("gamma_", length(records), ".rds")))
  write.csv(do.call(rbind, records), file.path(output, "checks.csv"), row.names = FALSE)
  print(record)
}
