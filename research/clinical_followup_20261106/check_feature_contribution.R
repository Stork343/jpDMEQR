library(data.table)
library(yaml)
source("research/clinical_search_20261105/panel.R")
folder <- "results/preflight/clinical_followup_20261106"
register <- fread(file.path(folder, "screening_register.csv"))[is.finite(cq_loss)]
result <- list()
for (i in seq_len(nrow(register))) {
  study <- read_yaml(register$study_settings[i])
  method <- read_yaml(register$method_settings[i])
  root <- file.path(method$output_root, study$study_id)
  panel <- cp_panel(study$input_directory)
  records <- list()
  active <- numeric(study$outer_folds)
  for (fold in seq_len(study$outer_folds)) {
    rows <- fread(file.path(study$output_root, paste0("outer_", fold), "evaluation.csv"))
    stopifnot(all(rows$role == "discovery_evaluation"))
    quantiles <- as.numeric(unlist(method$quantiles))
    full <- clinical <- matrix(NA_real_, nrow(rows), length(quantiles))
    for (j in seq_along(quantiles)) {
      tau <- quantiles[j]
      path <- file.path(root, paste0("outer_",fold), sprintf("tau%03d", round(tau*100)))
      fit <- readRDS(file.path(path, "fit.rds"))
      stopifnot(!length(intersect(fit$recipe$training_strains, rows$strain_id)))
      X <- cp_features(rows, panel, fit$recipe)$X[, names(fit$beta), drop = FALSE]
      base <- fit$clinical_columns
      stopifnot(all(base %in% colnames(X)), "(Intercept)" %in% base)
      full[, j] <- drop(X %*% fit$beta) * fit$recipe$y_scale + fit$recipe$y_center
      clinical[, j] <- drop(X[, base, drop = FALSE] %*% fit$beta[base]) * fit$recipe$y_scale + fit$recipe$y_center
      stored <- fread(file.path(path, "prediction.csv"))
      stopifnot(identical(stored$sample_id, rows$sample_id), max(abs(stored$prediction-full[, j])) < 1e-8)
      if (tau == study$decision_quantile) active[fold] <- sum(abs(fit$beta[!names(fit$beta) %in% base]) > 1e-7)
    }
    j <- match(study$decision_quantile,quantiles)
    a <- pmax(0, apply(full, 1L, sort)[j, ])
    b <- pmax(0, apply(clinical, 1L, sort)[j, ])
    loss <- function(prediction) {
      residual <- rows$y-prediction
      residual*(study$decision_quantile-as.numeric(residual<0))
    }
    records[[fold]] <- data.table(patient = rows$strain_id, full = loss(a), fixed_clinical = loss(b),
                                  absolute_prediction_change = abs(a-b))
  }
  patient <- rbindlist(records)[, lapply(.SD, mean), by = patient]
  full_loss <- mean(patient$full)
  ablated_loss <- mean(patient$fixed_clinical)
  stopifnot(nrow(patient) == register$patients[i], abs(full_loss-register$cq_loss[i]) < 1e-8)
  result[[i]] <- data.table(study_id = study$study_id, patients = nrow(patient), tau = study$decision_quantile,
    nonzero_extra_min = min(active), nonzero_extra_median = median(active), nonzero_extra_max = max(active),
    full_loss, fixed_clinical_loss = ablated_loss, gain_from_extra_terms_at_fixed_coefficients = 1-full_loss/ablated_loss,
    mean_absolute_prediction_change = mean(patient$absolute_prediction_change), reserved_outcomes_used = FALSE)
}
fwrite(rbindlist(result), file.path(folder,"discovery_feature_contribution.csv"))
