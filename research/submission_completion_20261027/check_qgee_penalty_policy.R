.libPaths(c(normalizePath("results/scratch/profile_pair_package_20261022/library"), .libPaths()))
source("research/submission_completion_20261027/application.R")
namespace <- asNamespace("geeVerse")
reference_kernel <- get(".fit_single_lambda", namespace)
reference_update <- get(".update_beta", namespace)
checked <- new.env(parent = namespace)
checked$.update_beta <- function(x_active, y, mu, beta_current, R, ...) {
  if (any(!is.finite(c(x_active, y, mu, beta_current, R)))) stop("Nonfinite QGEE state.")
  reference_update(x_active, y, mu, beta_current, R, ...)
}
environment(reference_kernel) <- checked
path <- "results/raw/application_comparison/submission_20261027_v2/job_0003/inner_fold_4_candidate_2"
input <- readRDS(file.path(path, "input.rds"))
output <- "results/scratch/submission_completion_20261027/qgee_penalty_policy_check"
if (dir.exists(output)) stop("The penalty-path policy check already exists.")
dir.create(output, recursive = TRUE)
reference <- submission_qgee(input$metadata, input$design, input$training, input$tau, input$seed,
  kernel = reference_kernel, allow_failed_penalties = TRUE)
candidate <- submission_qgee(input$metadata, input$design, input$training, input$tau, input$seed,
  allow_failed_penalties = TRUE)
saveRDS(list(reference = reference, candidate = candidate), file.path(output, "fits.rds"))
difference <- max(abs(reference$beta_prediction - candidate$beta_prediction))
if (difference > 1e-5 || reference$converged != candidate$converged ||
    abs(min(reference$lambda_path$hbic) - min(candidate$lambda_path$hbic)) > 1e-7)
  stop("The individual-penalty HBIC procedure differs from its maintained-kernel calculation.")
record <- data.frame(maximum_coefficient_difference = difference, converged = candidate$converged,
  failed_reference_penalties = sum(!is.finite(reference$lambda_path$hbic)),
  failed_candidate_penalties = sum(!is.finite(candidate$lambda_path$hbic)),
  selected_penalty = candidate$selected_penalty)
write.csv(record, file.path(output, "checks.csv"), row.names = FALSE)
print(record)
