library(data.table)
library(jsonlite)
library(yaml)
library(matrixStats)
source("research/score_validity_20261026/load_joint.R")
source("research/charite_application_20261105/features.R")
config <- read_yaml("config/application_charite_hb_20261105.yml")
root <- config$output_root
report <- config$report_directory
release <- read_json(file.path(root, "validation_release.json"), simplifyVector = TRUE)
rows <- fread(file.path(config$input_directory, "validation.csv"))
panel <- ch_panel(config$input_directory)
targets <- paste0("P_", release$target_proteins)
nt <- length(targets)
records <- list()
for (multiplier in unlist(config$inference$bandwidth_multipliers)) {
  saved <- readRDS(file.path(root, "inference", sprintf("bandwidth_%03d.rds", as.integer(multiplier * 100))))
  features <- ch_features(rows, panel, saved$recipe)
  features$X <- features$X[, saved$columns, drop = FALSE]
  data <- sb_data(features$y, features$X, features$Z, features$cluster_id)
  conversion <- saved$recipe$y_scale / saved$recipe$scale[targets]
  reconstructed <- estimate <- list()
  for (value in unlist(config$quantiles)) {
    fit <- saved$fits[[as.character(value)]]
    scores <- list()
    for (part in c("full", "half")) {
      fitted <- fit[[part]]
      h <- fit$diagnostics$h / if (part == "half") 2 else 1
      profile <- ds_profile(data, fitted$beta, value, h,
        matrix(config$inference$working_lambda), fitted$gamma)
      score <- t(vapply(data$rows, function(index) {
        -colSums(features$X[index, , drop = FALSE] *
          as.numeric(profile$psi[index]) * data$weight[index])
      }, numeric(ncol(features$X))))
      stopifnot(identical(rownames(score), saved$cluster_ids))
      precision <- solve(fitted$hessian)[, seq_len(nt), drop = FALSE]
      scores[[part]] <- score %*% precision
    }
    independent <- sweep(scores$full - 2 * scores$half, 2L, conversion, "*")
    coefficient <- (2 * fit$half$beta[seq_len(nt)] - fit$full$beta[seq_len(nt)]) * conversion
    reconstructed[[length(reconstructed) + 1L]] <- independent
    estimate[[length(estimate) + 1L]] <- coefficient
  }
  influence <- do.call(cbind, reconstructed)
  estimates <- unlist(estimate, use.names = FALSE)
  influence_error <- max(abs(influence - saved$influence))
  coefficient_error <- max(abs(estimates - saved$estimates))
  contrast_values <- drop(saved$contrasts %*% estimates)
  projected <- influence %*% t(saved$contrasts)
  standard_errors <- apply(projected, 2L, function(x) sqrt(sum((x - mean(x))^2)) / data$n)
  se_error <- max(abs(standard_errors - saved$labels$se))
  contrast_error <- max(abs(contrast_values - saved$labels$estimate))
  stopifnot(influence_error < 1e-9, coefficient_error < 1e-10,
    se_error < 1e-10, contrast_error < 1e-10)
  records[[length(records) + 1L]] <- data.table(multiplier, influence_error,
    coefficient_error, contrast_error, se_error, passed = TRUE)
}
fwrite(rbindlist(records), file.path(report, "inference_reconstruction_checks.csv"))
runtime <- lapply(release$model_directories, function(directory) {
  request <- read_json(file.path(directory, "request.json"), simplifyVector = TRUE)
  status <- read_json(file.path(directory, "status.json"), simplifyVector = TRUE)
  stopifnot(status$status == "completed", status$stage == "training")
  data.table(method = request$method, run_id = basename(directory),
    runtime_sec = status$runtime_sec, status = status$status)
})
fwrite(rbindlist(runtime), file.path(report, "training_runtime.csv"))
writeLines(trimws(capture.output(sessionInfo()), which = "right"), file.path(report, "audit_session.txt"))
write_json(list(status = "passed", inference_reconstructions = length(records),
  completed_training_jobs = length(runtime), patients = uniqueN(rows$strain_id),
  code_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE)),
  file.path(report, "r_audit_status.json"), auto_unbox = TRUE, pretty = TRUE)
print(rbindlist(records))
