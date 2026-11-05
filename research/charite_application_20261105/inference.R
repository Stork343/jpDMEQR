library(data.table)
library(jsonlite)
library(yaml)
library(matrixStats)
source("research/score_validity_20261026/load_joint.R")
source("research/joint_profile_manuscript_20261026/fixed_inference.R")
source("research/charite_application_20261105/features.R")
config <- read_yaml("config/application_charite_hb_20261105.yml")
root <- config$output_root
release <- read_json(file.path(root, "validation_release.json"), simplifyVector = TRUE)
stopifnot(release$status == "frozen")
targets <- as.character(release$target_proteins)
target_columns <- paste0("P_", targets)
rows <- fread(file.path(config$input_directory, "validation.csv"))
panel <- ch_panel(config$input_directory)
recipe <- readRDS(file.path(root, "recipes.rds"))$training$direction_recipe$base_recipe
stopifnot(!length(intersect(rows$strain_id, recipe$training_strains)), all(rows$role == "validation"))
factor <- recipe$y_scale / recipe$scale[target_columns]
stopifnot(all(is.finite(factor)), all(factor > 0))
build_design <- function(current) {
  features <- ch_features(current, panel, recipe)
  stopifnot(all(target_columns %in% colnames(features$X)))
  controls <- colnames(features$X)[colnames(features$X) == "(Intercept)" |
    startsWith(colnames(features$X), "C_") | colnames(features$X) %in% paste0("M_", targets)]
  selected <- target_columns
  weighted <- features$X * sqrt(features$weight)
  stopifnot(qr(weighted[, selected, drop = FALSE], tol = 1e-10)$rank == length(selected))
  for (column in controls) {
    candidate <- c(selected, column)
    if (qr(weighted[, candidate, drop = FALSE], tol = 1e-10)$rank == length(candidate)) selected <- candidate
  }
  features$X <- features$X[, selected, drop = FALSE]
  features$dropped_controls <- setdiff(controls, selected)
  features
}
features <- build_design(rows)
n <- features$n
q <- unlist(config$quantiles)
nt <- length(targets)
family_size <- nt * length(config$inference$contrasts)
Lambda <- matrix(as.numeric(config$inference$working_lambda))
output <- file.path(root, "inference")
if (dir.exists(output)) stop("推断目录已经存在")
dir.create(output)
coefficient_records <- contrast_records <- diagnostic_records <- list()
reference_influence <- reference_estimates <- NULL
contrast_matrix <- matrix(0, family_size, nt * length(q))
labels <- list()
cursor <- 0L
for (j in seq_len(nt)) for (k in 2:3) {
  cursor <- cursor + 1L
  contrast_matrix[cursor, j] <- 1
  contrast_matrix[cursor, (k - 1L) * nt + j] <- -1
  labels[[cursor]] <- data.table(protein = targets[j], contrast = if (k == 2L) "tau025_minus_tau050" else "tau025_minus_tau075")
}
labels <- rbindlist(labels)
for (multiplier in unlist(config$inference$bandwidth_multipliers)) {
  fits <- list()
  estimates <- influences <- list()
  for (value in q) {
    cat("Inference", multiplier, value, "\n")
    flush.console()
    fit <- jp_fit_fixed(features$y, features$X, features$Z, features$cluster_id,
      value, coordinates = seq_len(nt), Lambda = Lambda, h = n^(-0.3) * multiplier)
    stopifnot(max(fit$diagnostics$full_kkt, fit$diagnostics$half_kkt) < 1e-8,
      max(fit$diagnostics$full_nuisance_gradient, fit$diagnostics$half_nuisance_gradient) < 1e-8)
    estimates[[length(estimates) + 1L]] <- as.numeric(fit$table$beta_hat) * factor
    influences[[length(influences) + 1L]] <- sweep(fit$influence, 2L, factor, "*")
    coefficient_records[[length(coefficient_records) + 1L]] <- data.table(protein = targets,
      tau = value, multiplier, estimate = fit$table$beta_hat * factor, se = fit$table$se * factor,
      lower95 = fit$table$lower * factor, upper95 = fit$table$upper * factor,
      patients = n, records = nrow(rows), method = "JOINT_PROFILE", unit = "g_per_dl_per_log2_abundance")
    diagnostic_records[[length(diagnostic_records) + 1L]] <- cbind(data.table(tau = value, multiplier), as.data.table(fit$diagnostics))
    fits[[as.character(value)]] <- fit
  }
  estimate <- unlist(estimates, use.names = FALSE)
  influence <- do.call(cbind, influences)
  centered <- scale(influence, center = TRUE, scale = FALSE)
  covariance <- crossprod(centered) / n^2
  table <- copy(labels)
  table[, `:=`(multiplier = multiplier, estimate = drop(contrast_matrix %*% estimate),
    se = sqrt(diag(contrast_matrix %*% covariance %*% t(contrast_matrix))))]
  direct_se <- sqrt(colSums((centered %*% t(contrast_matrix))^2)) / n
  stopifnot(max(abs(direct_se - table$se)) < 1e-10, all(is.finite(table$se)), all(table$se > 0))
  critical <- qnorm(1 - config$validation$inference_alpha / (2 * family_size))
  table[, `:=`(lower = estimate - critical * se, upper = estimate + critical * se,
    p_value = 2 * pnorm(-abs(estimate / se)), patients = n, family_size = family_size,
    formal_validation_used = TRUE, reference_bandwidth = multiplier == 1)]
  table[, p_bonferroni := pmin(1, p_value * family_size)]
  table[, contrast_pass := lower > 0 | upper < 0]
  contrast_records[[length(contrast_records) + 1L]] <- table
  saveRDS(list(fits = fits, estimates = estimate, influence = influence, covariance = covariance,
    contrasts = contrast_matrix, labels = table, recipe = recipe, columns = colnames(features$X),
    cluster_ids = names(sb_data(features$y, features$X, features$Z, features$cluster_id)$rows)),
    file.path(output, sprintf("bandwidth_%03d.rds", as.integer(multiplier * 100))))
  if (multiplier == 1) {
    reference_influence <- influence
    reference_estimates <- estimate
  }
}
gram <- crossprod(features$X, features$X * features$weight) / n
mean_beta <- drop(qr.solve(features$X * sqrt(features$weight), features$y * sqrt(features$weight)))
data <- sb_data(features$y, features$X, features$Z, features$cluster_id)
mean_scores <- rowsum(features$X * ((features$y - drop(features$X %*% mean_beta)) * features$weight), data$group, reorder = FALSE)
mean_influence <- mean_scores %*% solve(gram)[, seq_len(nt), drop = FALSE]
mean_influence <- sweep(mean_influence, 2L, factor, "*")
mean_se <- sqrt(colSums(scale(mean_influence, center = TRUE, scale = FALSE)^2)) / n
mean_estimate <- mean_beta[seq_len(nt)] * factor
mean_table <- data.table(protein = targets, estimate = mean_estimate, se = mean_se,
  lower95 = mean_estimate - qnorm(0.975) * mean_se, upper95 = mean_estimate + qnorm(0.975) * mean_se,
  patients = n, method = "weighted_mean", interval_role = "descriptive")
ids <- sort(unique(rows$strain_id))
deletions <- matrix(NA_real_, length(ids), length(reference_estimates), dimnames = list(ids, NULL))
for (i in seq_along(ids)) {
  current <- build_design(rows[strain_id != ids[i]])
  for (k in seq_along(q)) {
    fit <- jp_fit_fixed(current$y, current$X, current$Z, current$cluster_id, q[k],
      coordinates = seq_len(nt), Lambda = Lambda, h = n^(-0.3))
    deletions[i, (k - 1L) * nt + seq_len(nt)] <- fit$table$beta_hat * factor
  }
}
stopifnot(all(is.finite(deletions)))
jackknife_covariance <- (n - 1) / n * crossprod(scale(deletions, center = TRUE, scale = FALSE))
jackknife <- copy(labels)
jackknife[, `:=`(estimate = drop(contrast_matrix %*% reference_estimates),
  se_jackknife = sqrt(diag(contrast_matrix %*% jackknife_covariance %*% t(contrast_matrix))))]
critical <- qnorm(1 - config$validation$inference_alpha / (2 * family_size))
jackknife[, `:=`(lower_jackknife = estimate - critical * se_jackknife,
  upper_jackknife = estimate + critical * se_jackknife)]
saveRDS(list(deletions = deletions, covariance = jackknife_covariance), file.path(output, "patient_deletions.rds"))
fwrite(rbindlist(coefficient_records), file.path(config$report_directory, "profile_coefficients.csv"))
fwrite(rbindlist(contrast_records), file.path(config$report_directory, "profile_contrasts.csv"))
fwrite(rbindlist(diagnostic_records), file.path(config$report_directory, "inference_diagnostics.csv"))
fwrite(mean_table, file.path(config$report_directory, "mean_associations.csv"))
fwrite(jackknife, file.path(config$report_directory, "jackknife_contrasts.csv"))
write_json(list(status = "completed", patients = n, records = nrow(rows), dimension = ncol(features$X),
  target_proteins = targets, clinical_columns = setdiff(colnames(features$X), target_columns),
  dropped_unidentified_controls = features$dropped_controls, family_size = family_size,
  inference_alpha = config$validation$inference_alpha, formal_validation_used = TRUE,
  manuscript_adoption = FALSE), file.path(config$report_directory, "inference_status.json"),
  auto_unbox = TRUE, pretty = TRUE, digits = NA)
print(rbindlist(contrast_records)[multiplier == 1])
