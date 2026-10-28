.libPaths(c(normalizePath("results/scratch/profile_pair_package_20261022/library"), .libPaths()))
source("research/submission_completion_20261027/application.R")
root <- "results/raw/application_comparison/submission_20261027_v2/job_0003"
output <- "results/scratch/submission_completion_20261027/qgee_lapack_path_checks"
if (dir.exists(output)) stop("The QGEE SVD comparison already exists.")
dir.create(output, recursive = TRUE)
records <- list()
for (j in 1:3) {
  path <- file.path(root, paste0("inner_fold_1_candidate_", j))
  input <- readRDS(file.path(path, "input.rds"))
  reference <- readRDS(file.path(path, "result.rds"))
  started <- proc.time()[[3L]]
  candidate <- submission_qgee(input$metadata, input$design, input$training, input$tau, input$seed,
    kernel = submission_qgee_kernel(submission_qgee_inverse_lapack))
  elapsed <- proc.time()[[3L]] - started
  difference <- max(abs(reference$beta_prediction - candidate$beta_prediction))
  objective_difference <- abs(min(reference$lambda_path$hbic) - min(candidate$lambda_path$hbic))
  saveRDS(candidate, file.path(output, paste0("candidate_", j, ".rds")))
  record <- data.frame(candidate = j, p = length(candidate$beta_prediction), elapsed = elapsed,
    coefficient_difference = difference, objective_difference = objective_difference,
    reference_penalty = reference$selected_penalty, selected_penalty = candidate$selected_penalty,
    converged = candidate$converged)
  records[[j]] <- record
  write.csv(do.call(rbind, records), file.path(output, "checks.csv"), row.names = FALSE)
  if (reference$converged != candidate$converged || difference > 1e-5 || objective_difference > 1e-7)
    stop("The QGEE application fit changed across SVD algorithms.")
  print(record)
}
