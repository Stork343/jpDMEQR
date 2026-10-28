.libPaths(c(normalizePath("results/scratch/profile_pair_package_20261022/library"), .libPaths()))
source("research/submission_completion_20261027/application.R")
source("research/submission_completion_20261027/qgee_lapack.R")
namespace <- asNamespace("geeVerse")
captured <- readRDS("results/scratch/submission_completion_20261027/qgee_inverse_diagnostic/last_matrix.rds")
reference_action <- get("geninv", namespace)(captured$matrix) %*% captured$rhs
candidate_action <- submission_qgee_inverse_lapack(captured$matrix, captured$rhs)
difference <- max(abs(candidate_action - reference_action))
if (difference > 1e-7) stop("The captured QGEE inverse action differs from its package definition.")
output <- "results/scratch/submission_completion_20261027/qgee_lapack_check"
if (dir.exists(output)) stop("The LAPACK reference check already exists.")
dir.create(output, recursive = TRUE)
write.csv(data.frame(dimension = nrow(captured$matrix), difference = difference),
  file.path(output, "captured_matrix_check.csv"), row.names = FALSE)
reference_kernel <- get(".fit_single_lambda", namespace)
reference_update <- get(".update_beta", namespace)
checked <- new.env(parent = namespace)
checked$.update_beta <- function(x_active, y, mu, beta_current, R, ...) {
  if (any(!is.finite(c(x_active, y, mu, beta_current, R)))) stop("Nonfinite QGEE state.")
  reference_update(x_active, y, mu, beta_current, R, ...)
}
environment(reference_kernel) <- checked
input <- readRDS("results/raw/application_comparison/submission_20261027_v2/job_0267/inner_fold_3_candidate_2/input.rds")
started <- proc.time()[[3L]]
reference <- submission_qgee(input$metadata, input$design, input$training, input$tau, input$seed,
  kernel = reference_kernel, allow_failed_penalties = TRUE)
reference_seconds <- proc.time()[[3L]] - started
started <- proc.time()[[3L]]
candidate <- submission_qgee(input$metadata, input$design, input$training, input$tau, input$seed,
  kernel = submission_qgee_kernel(submission_qgee_inverse_lapack), allow_failed_penalties = TRUE)
candidate_seconds <- proc.time()[[3L]] - started
coefficient_difference <- max(abs(reference$beta_prediction - candidate$beta_prediction))
if (coefficient_difference > 1e-5 || reference$converged != candidate$converged ||
    abs(min(reference$lambda_path$hbic) - min(candidate$lambda_path$hbic)) > 1e-7)
  stop("The complete LAPACK penalty path differs from its maintained-kernel calculation.")
saveRDS(list(reference = reference, candidate = candidate), file.path(output, "fits.rds"))
record <- data.frame(maximum_coefficient_difference = coefficient_difference,
  converged = candidate$converged, reference_seconds = reference_seconds, candidate_seconds = candidate_seconds,
  failed_reference_penalties = sum(!is.finite(reference$lambda_path$hbic)),
  failed_candidate_penalties = sum(!is.finite(candidate$lambda_path$hbic)))
write.csv(record, file.path(output, "checks.csv"), row.names = FALSE)
print(record)
