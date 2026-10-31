library(data.table)
library(yaml)
config <- read_yaml("config/application_soynam_20261031.yml")
arguments <- commandArgs(trailingOnly = TRUE)
mode <- if (length(arguments)) arguments[1L] else "all"
if (!mode %in% c("all", "reference")) stop("Unknown inference summary mode.")
registry <- fread(file.path(config$results_directory, "registry.csv"))
if (mode == "reference") registry <- registry[setting == "reference"]
root <- file.path(config$results_directory, "inference", config$execution_version)
destination <- if (mode == "reference") "results/preflight/soynam_20261031/reference" else "results/tables/soynam_20261031"
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
complete <- file.exists(file.path(root, registry$case_id, "completion.json"))
if (!all(complete)) stop("Incomplete registered fits: ", paste(registry$case_id[!complete], collapse = ", "))
estimates <- rbindlist(lapply(registry$case_id, function(id) fread(file.path(root, id, "inference.csv"))))
diagnostics <- rbindlist(lapply(registry$case_id, function(id) fread(file.path(root, id, "diagnostics.csv"), integer64 = "double")))
influences <- list()
betas <- standard_errors <- numeric()
labels <- character()
ids <- NULL
for (tau in unlist(config$quantiles)) {
  fit <- readRDS(file.path(root, sprintf("reference_q%02d", round(100 * tau)), "fit.rds"))
  current_ids <- unlist(lapply(fit$prepared$folds, `[[`, "validation_subjects"), use.names = FALSE)
  if (anyDuplicated(current_ids) || length(current_ids) != fit$prepared$n) stop("Invalid influence identifiers.")
  if (is.null(ids)) ids <- sort(current_ids)
  if (!setequal(ids, current_ids)) stop("The reference quantiles use different lines.")
  for (k in seq_along(config$coordinates)) {
    influence <- fit$details[[k]]$influence[match(ids, current_ids)] * config$response_divisor
    reconstructed <- sqrt(sum(influence^2)) / length(ids)
    reported <- fit$table$se[k] * config$response_divisor
    if (abs(reconstructed - reported) > 1e-10) stop("The saved influence does not reproduce its standard error.")
    influences[[length(influences) + 1L]] <- influence
    betas <- c(betas, fit$table$beta_tilde[k] * config$response_divisor)
    standard_errors <- c(standard_errors, reported)
    labels <- c(labels, paste0(config$coordinates[[k]], "_q", sprintf("%02d", round(100 * tau))))
  }
  rm(fit)
  gc(FALSE)
}
V <- do.call(cbind, influences)
covariance <- crossprod(V) / length(ids)^2
dimnames(covariance) <- list(labels, labels)
if (max(abs(sqrt(diag(covariance)) - standard_errors)) > 1e-10) stop("Covariance diagonal mismatch.")
C <- diag(6)
rownames(C) <- labels
for (year in 1:2) {
  contrast <- numeric(6)
  contrast[c(year, year + 4L)] <- c(-1, 1)
  C <- rbind(C, contrast)
  rownames(C)[nrow(C)] <- paste0("upper_minus_lower_", config$years[[year]])
}
for (q in 1:3) {
  contrast <- numeric(6)
  contrast[c(2L * q - 1L, 2L * q)] <- c(-1, 1)
  C <- rbind(C, contrast)
  rownames(C)[nrow(C)] <- paste0("year_2014_minus_2013_q", sprintf("%02d", 100 * config$quantiles[[q]]))
}
if (nrow(C) != config$family_size) stop("Inference family size mismatch.")
theta <- drop(C %*% betas)
se <- sqrt(diag(C %*% covariance %*% t(C)))
critical <- qnorm(1 - (1 - config$confidence) / (2 * config$family_size))
contrasts <- data.table(parameter = rownames(C), estimate = theta, se = se,
  lower = theta - critical * se, upper = theta + critical * se,
  p_unadjusted = 2 * pnorm(-abs(theta / se)),
  p_bonferroni = pmin(1, config$family_size * 2 * pnorm(-abs(theta / se))))
fwrite(estimates, file.path(destination, "inference.csv"))
fwrite(diagnostics, file.path(destination, "diagnostics.csv"))
fwrite(contrasts, file.path(destination, "reference_contrasts.csv"))
fwrite(data.table(parameter = labels, covariance), file.path(destination, "reference_covariance.csv"))
saveRDS(list(cluster_id = ids, influences = V, covariance = covariance, contrasts = C),
  file.path(config$results_directory, "reference_joint_inference.rds"))
print(contrasts)
