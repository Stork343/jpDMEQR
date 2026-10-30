source("research/submission_completion_20261027/load_comparators.R")
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(2L, 3L)) stop("Provide an input directory, candidate and optional solver.")
root <- args[1L]
candidate <- as.integer(args[2L])
algorithm <- if (length(args) == 3L) args[3L] else "qp_certified"
if (algorithm == "constraint_generation") source("research/submission_completion_20261027/clime_constraint_generation.R")
if (algorithm == "compiled_certified") source("research/submission_completion_20261027/load_compiled_clime.R")
input <- readRDS(file.path(root, "inputs.rds"))
gamma <- input$cells$gamma[candidate]
output <- file.path(root, sprintf("%s_matrix_%02d.rds", algorithm, candidate))
globalCallingHandlers(error = function(condition) saveRDS(condition, paste0(output, ".failure")))
started <- proc.time()[[3L]]
computed <- solve_clime_matrix_iid_v2(input$H, gamma, solver_opts = list(algorithm = algorithm))
elapsed <- proc.time()[[3L]] - started
if (!computed$feasible) stop("The previously feasible matrix was not recovered.")
reference <- input$precision_records[[candidate]]$result
coordinates <- match(rownames(reference$omega_raw_targets), rownames(computed$omega_raw))
errors <- c(raw = max(abs(computed$omega_raw[coordinates, ] - reference$omega_raw_targets)),
  symmetric = max(abs(computed$omega_symmetric[coordinates, ] - reference$omega_symmetric_targets)),
  objective = max(abs(rowSums(abs(computed$omega_raw)) - reference$raw_row_l1)))
loss <- sum((profile_dense_crossprod_v2(input$validation_H, computed$omega_symmetric) - diag(nrow(input$H)))^2) / nrow(input$H)
errors <- c(errors, loss = abs(loss - input$cells$loss[candidate]))
saveRDS(list(result = computed, elapsed = elapsed, errors = errors, loss = loss), output)
print(list(gamma = gamma, elapsed = elapsed, errors = errors, loss = loss))
if (max(errors) > 2e-5) stop("The complete matrix or validation loss changed.")
