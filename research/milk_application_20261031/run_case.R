library(data.table)
library(yaml)
config <- read_yaml("config/application_milk_20261031.yml")
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 1L) stop("Provide exactly one registered case identifier.")
registry <- fread(file.path(config$results_directory, "registry.csv"))
job <- registry[case_id == arguments[1L]]
if (nrow(job) != 1L) stop("The case identifier is not unique in the registry.")
directory <- file.path(config$results_directory, config$execution_version, job$case_id)
if (dir.exists(directory)) stop("The case directory already exists.")
dir.create(directory, recursive = TRUE)
backup <- file.path(config$backup_directory, config$execution_version, job$case_id)
dir.create(backup, recursive = TRUE, showWarnings = FALSE)
started <- Sys.time()
options(error = function() {
  failure <- list(case_id = job$case_id, message = geterrmessage(), calls = lapply(sys.calls(), deparse),
    started = as.character(started), failed = as.character(Sys.time()))
  saveRDS(failure, file.path(directory, "failure.rds"))
  file.copy(file.path(directory, "failure.rds"), backup, overwrite = FALSE)
  quit(save = "no", status = 1L, runLast = FALSE)
})
source("research/score_validity_20261026/load_joint.R")
source("research/milk_application_20261031/nodewise.R")
source("research/milk_application_20261031/weighted_solver.R")
ds_fit_small <- milk_fit_small
milk_design_cache <- readRDS(file.path(config$results_directory, "design_cache.rds"))
pc_nodewise <- milk_nodewise_cached
sb_lasso_supports <- milk_curvature_supports
cohort <- readRDS(file.path(config$results_directory, "cohorts.rds"))[[job$cohort]]
request <- list(job = job, configuration = config,
  source_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE), started = as.character(started))
saveRDS(request, file.path(directory, "request.rds"))
cat(format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE), job$case_id, "started\n")
flush.console()
fit <- ds_fit_profile(cohort$y, cohort$X, cohort$Z, cohort$cluster_id,
  tau = job$tau, coordinates = config$coordinate, Lambda = matrix(job$lambda, 1L, 1L),
  penalty_factor = cohort$penalty_factor, h = job$h, confidence = config$confidence,
  keep_prepared = TRUE, budget_multiplier = config$budget_multiplier)
if (fit$method_version != config$method_version) stop("The fitted method version differs from the registry.")
table <- as.data.table(fit$table)
critical <- qnorm(1 - (1 - config$confidence) / (2 * config$family_size))
table[, `:=`(case_id = job$case_id, setting = job$setting, cohort = job$cohort,
  tau = job$tau, lambda = job$lambda, h = job$h, n = job$n, p = job$p,
  simultaneous_lower = beta_tilde - critical * se,
  simultaneous_upper = beta_tilde + critical * se, runtime_seconds = fit$timing[["total"]])]
diagnostics <- rbindlist(lapply(1:2, function(fold) {
  detail <- fit$details[[1L]]$fits[[fold]]
  model <- fit$prepared$folds[[fold]]$models[[1L]]
  rbindlist(lapply(c("full", "half"), function(bandwidth) {
    current <- detail[[bandwidth]]
    fitted_bandwidth <- if (bandwidth == "full") job$h else job$h / 2
    profile <- ds_profile(model$data, current$beta, job$tau, fitted_bandwidth,
      matrix(job$lambda, 1L, 1L), current$gamma)
    eigenvalues <- eigen(current$hessian, symmetric = TRUE, only.values = TRUE)$values
    basis_eigenvalues <- eigen(current$basis_hessian, symmetric = TRUE, only.values = TRUE)$values
    data.table(case_id = job$case_id, fold = fold, bandwidth = bandwidth,
      selected_dimension = length(model$support), estimation_cows = model$data$n,
      kkt = current$kkt, iterations = current$iterations,
      nuisance_gradient = max(profile$gradient_max),
      hessian_minimum = min(eigenvalues), hessian_maximum = max(eigenvalues),
      hessian_condition = max(eigenvalues) / min(eigenvalues),
      basis_condition = max(basis_eigenvalues) / min(basis_eigenvalues),
      prediction_identity_error = current$prediction_error,
      cluster_score_variance = var(current$projected_scores))
  }))
}))
if (any(!is.finite(diagnostics$hessian_condition)) || any(diagnostics$hessian_minimum <= 0) ||
    any(diagnostics$kkt > 1e-8) || any(diagnostics$nuisance_gradient > 1e-8))
  stop("The fit failed its numerical checks.")
saveRDS(fit, file.path(directory, "fit.rds"), compress = FALSE)
fwrite(table, file.path(directory, "inference.csv"))
fwrite(diagnostics, file.path(directory, "diagnostics.csv"))
completion <- list(case_id = job$case_id, status = "complete", started = as.character(started),
  completed = as.character(Sys.time()), wall_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")))
jsonlite::write_json(completion, file.path(directory, "completion.json"), auto_unbox = TRUE, pretty = TRUE)
if (!all(file.copy(list.files(directory, full.names = TRUE), backup, overwrite = FALSE))) stop("Result backup failed.")
print(table)
cat(format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE), job$case_id, "completed\n")
