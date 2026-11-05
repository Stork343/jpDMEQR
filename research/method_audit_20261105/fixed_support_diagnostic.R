library(data.table)
library(jsonlite)
library(yaml)
library(matrixStats)
source("research/score_validity_20261026/load_joint.R")
source("research/joint_profile_manuscript_20261026/conic_search.R")
source("research/clinical_search_20261105/panel.R")
source("research/empirical_selection_20261101/prediction_core.R")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
config <- read_yaml(args[1L])
tau <- config$decision_quantile
panel <- cp_panel(config$input_directory)
root <- file.path("results/raw/method_audit_20261105", config$study_id)
if (dir.exists(root)) stop("诊断目录已经存在")
dir.create(root, recursive = TRUE)
report <- file.path("results/preflight/method_audit_20261105", config$study_id)
dir.create(report, recursive = TRUE)
records <- diagnostics <- list()

fit_basis <- function(features, columns, penalty, nuisance) {
  X <- features$X[, columns, drop = FALSE]
  weight <- sqrt(features$weight / features$n)
  decomposition <- qr(X * weight, tol = 1e-10)
  keep <- decomposition$pivot[seq_len(decomposition$rank)]
  selected <- columns[keep]
  X <- features$X[, selected, drop = FALSE]
  decomposition <- qr(X * weight, tol = 1e-10)
  stopifnot(decomposition$rank == ncol(X))
  basis <- qr.Q(decomposition)
  triangular <- qr.R(decomposition)
  reconstruction <- qr.fitted(decomposition, features$X[, columns, drop = FALSE] * weight)
  span_error <- max(abs(reconstruction - features$X[, columns, drop = FALSE] * weight))
  stopifnot(span_error < 1e-8)
  Z <- matrix(as.numeric(nuisance), nrow(X), 1L)
  group <- as.integer(factor(features$cluster_id)) - 1L
  solution <- jp_joint_qp(basis / weight, Z, features$y, group, features$weight,
    matrix(penalty), numeric(ncol(X)), tau, numeric(ncol(X)), matrix(0, features$n, 1L),
    numeric(nrow(X)), gap_tolerance = 1e-8)
  beta <- setNames(numeric(ncol(features$X)), colnames(features$X))
  beta[selected[decomposition$pivot]] <- drop(backsolve(triangular, solution$theta))
  if (!nuisance) stopifnot(max(abs(solution$gamma)) < 1e-8)
  list(beta = beta, solution = solution, columns = selected, span_error = span_error)
}

for (fold in seq_len(config$outer_folds)) {
  outer <- file.path(config$output_root, paste0("outer_", fold))
  train_rows <- fread(file.path(outer, "training.csv"))
  test_rows <- fread(file.path(outer, "evaluation.csv"))
  stopifnot(all(train_rows$role == "training"), all(test_rows$role == "discovery_evaluation"),
    !length(intersect(train_rows$strain_id, test_rows$strain_id)))
  directory <- file.path(outer, "models", sprintf("JOINT_DIRECT_tau%03d", as.integer(tau * 100)))
  original <- readRDS(file.path(directory, "fit.rds"))
  train <- cp_features(train_rows, panel, original$recipe)
  test <- cp_features(test_rows, panel, original$recipe)
  setting <- as.data.table(original$grid)[original$choice]
  support <- unique(c(which(train$penalty_factor == 0),
    vapply(original$final_diagnostics[[1L]], `[[`, integer(1), "selected_coordinate")))
  clinical <- which(colnames(train$X) == "(Intercept)" | startsWith(colnames(train$X), "C_") |
    startsWith(colnames(train$X), "M_C_"))
  full <- unique(c(clinical, support))
  settings <- data.table(variant = c("profile_selected", "pinball_selected", "profile_full_clinical", "pinball_full_clinical"),
    nuisance = c(TRUE, FALSE, TRUE, FALSE), full_clinical = c(FALSE, FALSE, TRUE, TRUE))
  for (j in seq_len(nrow(settings))) {
    variant <- settings$variant[j]
    chosen <- if (settings$full_clinical[j]) full else support
    cat(config$study_id, fold, variant, "\n")
    flush.console()
    fitted <- fit_basis(train, chosen, setting$Lambda, settings$nuisance[j])
    training_prediction <- drop(train$X %*% fitted$beta) * original$recipe$y_scale + original$recipe$y_center
    calibration <- ep_calibrate(train_rows, training_prediction, tau, "strain_id", "condition")
    prediction <- drop(test$X %*% fitted$beta) * original$recipe$y_scale + original$recipe$y_center + calibration[test_rows$condition]
    error <- if (variant == "profile_selected") max(abs(prediction - original$prediction[, 1L])) else NA_real_
    if (variant == "profile_selected") stopifnot(error < 1e-6)
    records[[length(records) + 1L]] <- data.table(study_id = config$study_id, fold, variant,
      sample_id = test_rows$sample_id, strain_id = test_rows$strain_id, tau, y = test_rows$y,
      prediction = prediction, loss = ep_loss(test_rows$y, prediction, tau))
    diagnostics[[length(diagnostics) + 1L]] <- data.table(study_id = config$study_id, fold, variant,
      Lambda = setting$Lambda, requested_columns = length(chosen), rank = length(fitted$columns),
      gap = fitted$solution$duality_gap, beta_kkt = fitted$solution$beta_kkt,
      gamma_kkt = fitted$solution$gamma_kkt, span_error = fitted$span_error, original_prediction_error = error)
    saveRDS(list(fit = fitted, calibration = calibration, original_support = support,
      original_directory = directory, tau = tau, variant = variant),
      file.path(root, sprintf("fold_%d_%s.rds", fold, variant)))
  }
}
all <- rbindlist(records)
patient <- all[, .(loss = mean(loss)), by = .(study_id, fold, variant, strain_id)]
summary <- patient[, .(patients = .N, loss = mean(loss), descriptive_se = sd(loss) / sqrt(.N)), by = .(study_id, variant)]
wide <- dcast(patient, study_id + strain_id ~ variant, value.var = "loss")
comparisons <- rbindlist(lapply(setdiff(unique(patient$variant), "profile_selected"), function(value) {
  difference <- wide$profile_selected - wide[[value]]
  data.table(study_id = config$study_id, variant = value, absolute_gain = mean(difference),
    relative_gain = mean(difference) / mean(wide$profile_selected), patients = nrow(wide))
}))
fwrite(all, file.path(root, "predictions.csv"))
fwrite(summary, file.path(report, "fixed_support_losses.csv"))
fwrite(comparisons, file.path(report, "fixed_support_comparisons.csv"))
fwrite(rbindlist(diagnostics), file.path(report, "solver_checks.csv"))
write_json(list(status = "completed", role = "exploratory_method_diagnostic", study_id = config$study_id,
  tau = tau, formal_validation_used = FALSE, manuscript_adoption = FALSE), file.path(report, "status.json"),
  auto_unbox = TRUE, pretty = TRUE)
print(summary)
print(comparisons)
