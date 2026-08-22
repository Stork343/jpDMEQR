#!/usr/bin/env Rscript

full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else
  "scripts/application/03_fit_GSE65391.R"
root <- if (file.exists("R/application_v2.R")) "." else
  normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", envir = environment())
source_v2_module(root, "metrics_v2", envir = environment())
source_v2_module(root, "application_v2", envir = environment())

parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!grepl("^--", arg)) next
    pair <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    out[[pair[1]]] <- if (length(pair) > 1L) paste(pair[-1], collapse = "=") else "true"
  }
  out
}
as_bool <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  tolower(as.character(x)) %in% c("true", "1", "yes", "y")
}
as_num_vec <- function(x, default) {
  if (is.null(x) || !nzchar(x)) return(default)
  as.numeric(strsplit(x, ",", fixed = TRUE)[[1]])
}
args <- parse_args(commandArgs(trailingOnly = TRUE))
mode <- args$mode %||% "smoke"
repetition <- as.integer(args$repetition %||% 1L)
fold_id <- as.integer(args$fold %||% 1L)
models <- toupper(strsplit(args$models %||% "A", ",", fixed = TRUE)[[1]])
taus <- as_num_vec(args$taus, if (mode == "smoke") 0.50 else c(0.25, 0.50, 0.75))
top_features <- as.integer(args$features %||% if (mode == "smoke") 100L else 5000L)
h_multiplier <- as.numeric(args$h_multiplier %||% 1)
lambda_beta_multiplier <- as.numeric(args$lambda_beta_multiplier %||% 1)
lambda_gamma <- as.numeric(args$lambda_gamma %||% 1)
random_slope <- as_bool(args$random_slope, FALSE)
run_inference <- as_bool(args$inference, FALSE)
max_iter <- as.integer(args$max_iter %||% if (mode == "smoke") 40L else 200L)
run_id <- args$run_id %||% paste0(
  format(Sys.time(), "%Y%m%dT%H%M%S"), "_rep", repetition, "_fold", fold_id, "_", mode
)

if (any(!is.finite(taus) | taus <= 0 | taus >= 1)) stop("Invalid --taus.")
if (!all(models %in% c("A", "B"))) stop("--models must contain only A and/or B.")

cohort_path <- file.path(root, "data", "derived", "GSE65391", "analysis_cohort.rds")
annotation_path <- file.path(root, "data", "interim", "GSE65391", "probe_annotation.csv")
if (!file.exists(cohort_path) || !file.exists(annotation_path)) {
  stop("Run 00_download, 01_build and 02_audit before fitting.")
}
obj <- readRDS(cohort_path)
if (!isTRUE(obj$audit$summary$all_gates_pass)) stop("OUTCOME_GATE_FAILED: fitting is blocked.")
metadata <- obj$eligible_metadata
expression_probe <- obj$expression[, rownames(metadata), drop = FALSE]
folds <- obj$subject_folds
fold_row <- folds$repetition == repetition & folds$fold == fold_id
if (!any(fold_row)) stop("Requested repetition/fold was not found.")
test_subjects <- as.character(folds$subject_id[fold_row])
train_rows <- !metadata$subject_id %in% test_subjects
test_rows <- !train_rows
if (!any(train_rows) || !any(test_rows)) stop("Empty training or test set.")

annotation <- read.csv(annotation_path, stringsAsFactors = FALSE, check.names = FALSE)
for (v in intersect(c("is_control", "is_multimapping"), names(annotation))) {
  annotation[[v]] <- as.logical(annotation[[v]])
  annotation[[v]][is.na(annotation[[v]])] <- FALSE
}
filtered <- remove_probe_rows_v2(expression_probe, annotation)
collapsed <- collapse_probes_highest_training_mad_v2(
  filtered$expression,
  filtered$annotation$gene_symbol,
  training_samples = rownames(metadata)[train_rows]
)
mad_filtered <- filter_top_training_mad_v2(
  collapsed$expression,
  training_samples = rownames(metadata)[train_rows],
  d = top_features
)
expression_gene <- mad_filtered$expression
X_gene_raw <- t(expression_gene[, rownames(metadata), drop = FALSE])
gene_scaled <- scale_from_training_v2(X_gene_raw, train_rows)
X_gene <- gene_scaled$X
colnames(X_gene) <- paste0("gene__", colnames(X_gene))

out_dir <- file.path(root, "results", "raw", "application", run_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(data.frame(
  sample_id = rownames(metadata), subject_id = metadata$subject_id,
  split = ifelse(train_rows, "train", "test"),
  stringsAsFactors = FALSE
), file.path(out_dir, "split.csv"), row.names = FALSE)
write.csv(data.frame(
  gene_symbol = sub("^gene__", "", colnames(X_gene)),
  training_mad = mad_filtered$mad[gene_scaled$keep],
  center = gene_scaled$center,
  scale = gene_scaled$scale,
  stringsAsFactors = FALSE
), file.path(out_dir, "gene_preprocessing.csv"), row.names = FALSE)

all_metrics <- list(); all_coef <- list(); all_pred <- list(); all_inf <- list()
metric_index <- coef_index <- pred_index <- inf_index <- 0L
for (model_id in models) {
  clinical <- prepare_clinical_design_split_v2(
    metadata, train_rows, include_sledai = identical(model_id, "B")
  )
  X_all <- cbind(clinical$X, X_gene)
  penalty_factor <- c(rep(0, ncol(clinical$X)), rep(1, ncol(X_gene)))
  Z_all <- if (random_slope) cbind(intercept = 1, centered_time = metadata$cumulative_time - clinical$time_center) else
    matrix(1, nrow(metadata), 1L, dimnames = list(NULL, "intercept"))

  for (tau in taus) {
    n_train_clusters <- length(unique(metadata$subject_id[train_rows]))
    h <- h_multiplier * n_train_clusters^(-3 / 10)
    lambda_beta <- lambda_beta_multiplier * sqrt(log(max(ncol(X_all), 2L)) / n_train_clusters)
    fit_start <- proc.time()[[3]]
    fit <- tryCatch(
      fit_profile_lasso_v2(
        y = metadata$log_c3[train_rows],
        X = X_all[train_rows, , drop = FALSE],
        Z = Z_all[train_rows, , drop = FALSE],
        cluster_id = metadata$subject_id[train_rows],
        tau = tau,
        h = h,
        lambda_beta = lambda_beta,
        lambda_gamma = lambda_gamma,
        penalty_factor = penalty_factor,
        control = list(
          max_iter = max_iter,
          beta_tol = 1e-6,
          kkt_tol = 1e-5,
          nuisance_control = list(reltol = 1e-11, maxit = 500L, grad_tol = 1e-8)
        )
      ),
      error = function(e) e
    )
    runtime <- proc.time()[[3]] - fit_start
    if (inherits(fit, "error")) {
      metric_index <- metric_index + 1L
      all_metrics[[metric_index]] <- data.frame(
        run_id = run_id, model = model_id, tau = tau, status = "failed",
        failure_message = conditionMessage(fit), n_train_subjects = n_train_clusters,
        n_train_visits = sum(train_rows), n_test_subjects = length(unique(metadata$subject_id[test_rows])),
        n_test_visits = sum(test_rows), p_total = ncol(X_all), p_genes = ncol(X_gene),
        pinball_loss = NA_real_, quantile_calibration = NA_real_, selected_genes = NA_integer_,
        runtime_sec = runtime, converged = FALSE, kkt_residual = NA_real_,
        profile_identity_error = NA_real_, max_nuisance_gradient = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }

    qhat <- as.numeric(X_all[test_rows, , drop = FALSE] %*% fit$beta)
    metric_index <- metric_index + 1L
    fit_status <- if (fit$converged) "ok" else "warning"
    all_metrics[[metric_index]] <- data.frame(
      run_id = run_id, model = model_id, tau = tau, status = fit_status,
      failure_message = fit$failure_stage %||% "", n_train_subjects = n_train_clusters,
      n_train_visits = sum(train_rows), n_test_subjects = length(unique(metadata$subject_id[test_rows])),
      n_test_visits = sum(test_rows), p_total = ncol(X_all), p_genes = ncol(X_gene),
      pinball_loss = pinball_loss_v2(metadata$log_c3[test_rows], qhat, tau),
      quantile_calibration = mean(metadata$log_c3[test_rows] <= qhat),
      selected_genes = sum(abs(fit$beta[grepl("^gene__", names(fit$beta))]) > 1e-8),
      runtime_sec = runtime, converged = fit$converged, kkt_residual = fit$kkt_residual,
      profile_identity_error = fit$components$profile_identity_error,
      max_nuisance_gradient = fit$components$max_nuisance_gradient,
      stringsAsFactors = FALSE
    )

    coef_index <- coef_index + 1L
    all_coef[[coef_index]] <- data.frame(
      run_id = run_id, model = model_id, tau = tau,
      coordinate = names(fit$beta), beta_hat = as.numeric(fit$beta),
      coordinate_type = ifelse(grepl("^gene__", names(fit$beta)), "gene", "clinical"),
      selected = abs(as.numeric(fit$beta)) > 1e-8,
      stringsAsFactors = FALSE
    )
    pred_index <- pred_index + 1L
    all_pred[[pred_index]] <- data.frame(
      run_id = run_id, model = model_id, tau = tau,
      sample_id = rownames(metadata)[test_rows], subject_id = metadata$subject_id[test_rows],
      observed_log_c3 = metadata$log_c3[test_rows], predicted_quantile = qhat,
      stringsAsFactors = FALSE
    )

    if (run_inference) {
      if (!requireNamespace("CVXR", quietly = TRUE)) stop("--inference=true requires CVXR.")
      # Ordinary pointwise inference is restricted to always-included
      # clinical coordinates. Gene-level inference requires the independent
      # split-sample workflow described in docs/EMPIRICAL_PROTOCOL_GSE65391.md.
      target_coordinates <- names(fit$beta)[!grepl("^gene__", names(fit$beta))]
      mu0 <- sqrt(log(max(ncol(X_all), 2L)) / (n_train_clusters * h)) + h^2
      inf <- debias_profile_coordinates_v2(
        fit, target_coordinates, mu_grid = c(0.5, 1, 2, 4) * mu0,
        ci_level = 0.95, df_correction = TRUE
      )
      inf_index <- inf_index + 1L
      all_inf[[inf_index]] <- cbind(
        data.frame(run_id = run_id, model = model_id, tau = tau, stringsAsFactors = FALSE),
        inf$table
      )
    }
  }
}

metrics <- do.call(rbind, all_metrics)
coefficients <- if (length(all_coef)) do.call(rbind, all_coef) else data.frame()
predictions <- if (length(all_pred)) do.call(rbind, all_pred) else data.frame()
inference <- if (length(all_inf)) do.call(rbind, all_inf) else data.frame()
write.csv(metrics, file.path(out_dir, "fold_metrics.csv"), row.names = FALSE)
if (nrow(coefficients)) write.csv(coefficients, file.path(out_dir, "coefficients.csv"), row.names = FALSE)
if (nrow(predictions)) write.csv(predictions, file.path(out_dir, "predictions.csv"), row.names = FALSE)
if (nrow(inference)) write.csv(inference, file.path(out_dir, "inference.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
saveRDS(list(metrics = metrics, coefficients = coefficients, predictions = predictions,
             inference = inference, clinical = clinical, gene_scaled = gene_scaled),
        file.path(out_dir, "analysis_objects.rds"), compress = "xz")
message("Application fold complete: ", out_dir)
print(metrics)
