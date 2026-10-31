library(data.table)
library(yaml)
config <- read_yaml("config/application_soynam_20261031.yml")
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 1L) stop("Provide one registered inference case.")
registry <- fread(file.path(config$results_directory, "registry.csv"))
job <- registry[case_id == arguments[1L]]
if (nrow(job) != 1L) stop("The requested inference case is not unique.")
directory <- file.path(config$results_directory, "inference", job$case_id)
if (dir.exists(directory)) stop("The requested inference directory already exists.")
dir.create(directory, recursive = TRUE)
backup <- file.path(config$backup_directory, "inference", job$case_id)
dir.create(backup, recursive = TRUE, showWarnings = FALSE)
started <- Sys.time()
options(error = function() {
  failure <- list(case_id = job$case_id, error = geterrmessage(), calls = lapply(sys.calls(), deparse),
    started = as.character(started), failed = as.character(Sys.time()))
  saveRDS(failure, file.path(directory, "failure.rds"))
  file.copy(file.path(directory, "failure.rds"), backup, overwrite = FALSE)
  quit(save = "no", status = 1L, runLast = FALSE)
})
source("research/soynam_application_20261031/engine.R")
input <- readRDS(file.path(config$results_directory, "input.rds"))
data <- soy_design(input, job$cohort)
index <- cluster_index_v2(data$cluster_id)$rows
assignment <- pc_folds(index, 2L)
cache <- file.path(config$results_directory, "cache")
dir.create(cache, showWarnings = FALSE)
soy_context <- list(case_id = job$case_id, cohort = job$cohort, directory = directory,
  cache = cache, selection_subjects = lapply(1:2, function(fold) names(index)[assignment != fold]))
request <- list(job = job, configuration = config, n = data$n, p = data$p,
  source_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  started = as.character(started), session = sessionInfo())
saveRDS(request, file.path(directory, "request.rds"))
soy_progress("complete_high_dimensional_fit", list(clusters = data$n, variables = data$p))
fit <- ds_fit_profile(data$y, data$X, data$Z, data$cluster_id, job$tau,
  coordinates = unlist(config$coordinates), Lambda = matrix(job$lambda, 1L, 1L),
  penalty_factor = data$penalty_factor, h = job$h_multiplier * data$n^(-0.3),
  confidence = config$confidence, keep_prepared = TRUE, budget_multiplier = config$budget_multiplier)
table <- as.data.table(fit$table)
critical <- qnorm(1 - (1 - config$confidence) / (2 * config$family_size))
table[, `:=`(case_id = job$case_id, setting = job$setting, tau = job$tau,
  lambda = job$lambda, h = fit$h, n = data$n, p = data$p,
  estimate_bu_ac = config$response_divisor * beta_tilde,
  se_bu_ac = config$response_divisor * se,
  simultaneous_lower = config$response_divisor * (beta_tilde - critical * se),
  simultaneous_upper = config$response_divisor * (beta_tilde + critical * se),
  runtime_seconds = fit$timing[["total"]])]
diagnostics <- rbindlist(lapply(seq_along(config$coordinates), function(coordinate) {
  rbindlist(lapply(1:2, function(fold) {
    model <- fit$prepared$folds[[fold]]$models[[coordinate]]
    detail <- fit$details[[coordinate]]$fits[[fold]]
    rbindlist(lapply(c("full", "half"), function(bandwidth) {
      current <- detail[[bandwidth]]
      h <- if (bandwidth == "full") fit$h else fit$h / 2
      profile <- ds_profile(model$data, current$beta, job$tau, h, matrix(job$lambda, 1L, 1L), current$gamma)
      eigenvalues <- eigen(current$basis_hessian, symmetric = TRUE, only.values = TRUE)$values
      data.table(case_id = job$case_id, coordinate = unlist(config$coordinates)[coordinate],
        fold = fold, bandwidth = bandwidth, dimension = length(model$support),
        kkt = current$kkt, nuisance_gradient = max(profile$gradient_max),
        basis_minimum = min(eigenvalues), basis_condition = max(eigenvalues) / min(eigenvalues),
        prediction_identity_error = current$prediction_error, iterations = current$iterations)
    }))
  }))
}))
if (any(diagnostics$kkt > 1e-8) || any(diagnostics$nuisance_gradient > 1e-8) ||
    any(diagnostics$basis_minimum <= 0) || any(!is.finite(table$se_bu_ac)))
  stop("The empirical fit failed its numerical checks.")
saveRDS(fit, file.path(directory, "fit.rds"), compress = FALSE)
fwrite(table, file.path(directory, "inference.csv"))
fwrite(diagnostics, file.path(directory, "diagnostics.csv"))
jsonlite::write_json(list(case_id = job$case_id, status = "complete", started = as.character(started),
  completed = as.character(Sys.time()), wall_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))),
  file.path(directory, "completion.json"), auto_unbox = TRUE, pretty = TRUE)
files <- list.files(directory, full.names = TRUE)
if (!all(file.copy(files, backup, overwrite = FALSE))) stop("The completed empirical result backup failed.")
print(table)
