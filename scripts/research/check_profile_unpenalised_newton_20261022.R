root <- normalizePath(".", mustWork = TRUE)
.libPaths(c(file.path(root, "results", "scratch", "profile_pair_package_20261022", "library"), .libPaths()))
library(jpDMEQR)
Sys.setenv(JPDMEQR_CLUSTER_SOLVER = "compiled")
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", environment())
source_v2_module(root, "application_v2", environment())
source_v2_module(root, "simulation_v2", environment())
args <- commandArgs(trailingOnly = TRUE)
run_id <- if (length(args)) args[1] else "unpenalised_newton_20261022_v1"
output <- file.path(root, "results", "scratch", run_id)
if (dir.exists(output)) stop("The numerical validation directory already exists.")
dir.create(output, recursive = TRUE)
population_dir <- file.path(root, "data", "derived", "GSE65391", "execution_20261022_v1")
results <- list()
control <- list(max_iter = 1000L, max_backtrack = 50L, beta_tol = 1e-8,
                kkt_normalized_tol = 1e-7,
                nuisance_control = list(reltol = 1e-11, maxit = 400L, grad_tol = 1e-8))
for (model in c("A", "B")) {
  for (q in 1:2) {
    pop <- readRDS(file.path(population_dir, sprintf("matched_population_model_%s_q%d.rds", model, q)))
    y <- log(pop$metadata$c3)
    for (tau in c(.25, .5, .75)) {
      label <- sprintf("observed_%s_q%d_tau%02d", model, q, 100 * tau)
      run <- function(optimizer, beta_init = NULL) {
        started <- proc.time()[[3]]
        fit <- fit_profile_lasso_v2(
          y, pop$X, pop$Z, pop$metadata$subject_id, tau, 129^(-3 / 10), 0,
          lambda_gamma = diag(q), penalty_factor = rep(0, ncol(pop$X)),
          beta_init = beta_init,
          control = modifyList(control, list(unpenalised_optimizer = optimizer))
        )
        list(fit = fit, runtime_sec = proc.time()[[3]] - started)
      }
      newton <- run("newton_armijo")
      if (!isTRUE(newton$fit$converged)) stop("The complete observed design Newton fit did not converge: ", label)
      reference <- run("proximal_gradient", newton$fit$beta)
      cold <- run("proximal_gradient")
      saveRDS(list(newton = newton, reference_warm = reference, reference_zero = cold),
              file.path(output, paste0(label, ".rds")))
      row <- data.frame(
        case = label, n_clusters = 129L, n_visits = nrow(pop$X), p = ncol(pop$X),
        tau = tau, q = q, newton_converged = newton$fit$converged,
        reference_warm_converged = reference$fit$converged,
        reference_zero_converged = cold$fit$converged,
        newton_runtime_sec = newton$runtime_sec,
        reference_warm_runtime_sec = reference$runtime_sec,
        reference_zero_runtime_sec = cold$runtime_sec,
        newton_iterations = newton$fit$iterations,
        reference_zero_iterations = cold$fit$iterations,
        newton_kkt = newton$fit$kkt_normalized,
        reference_warm_kkt = reference$fit$kkt_normalized,
        reference_zero_kkt = cold$fit$kkt_normalized,
        max_coefficient_difference = max(abs(newton$fit$beta - reference$fit$beta)),
        objective_difference = abs(newton$fit$profile_objective - reference$fit$profile_objective),
        identity_error = newton$fit$components$profile_identity_error,
        scope = "numerical_computation_validation", stringsAsFactors = FALSE
      )
      if (!reference$fit$converged || row$max_coefficient_difference > 1e-6 ||
          row$objective_difference > 1e-10 || row$newton_kkt > 1e-7) {
        stop("A complete observed design numerical comparison failed: ", label)
      }
      results[[length(results) + 1L]] <- row
      write.csv(do.call(rbind, results), file.path(output, "observed_design_comparison.csv"), row.names = FALSE)
      print(row)
    }
  }
}
writeLines(capture.output(sessionInfo()), file.path(output, "sessionInfo.txt"))
