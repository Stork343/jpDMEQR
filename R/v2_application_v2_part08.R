prepare_gse65391_profile_workflow_stage_v2 <- function(metadata, expression_probe, annotation,
                                                       module_signatures, module_score_members,
                                                       module_rule, include_sledai,
                                                       coordinate_subset = NULL,
                                                       expected_coordinates = NULL,
                                                       random_slope = FALSE) {
  design <- prepare_gse65391_lowdim_design_v2(
    metadata, expression_probe, annotation, module_signatures,
    rep(TRUE, nrow(metadata)), module_rule, include_sledai = include_sledai,
    module_score_members = module_score_members
  )
  coordinates <- if (is.null(coordinate_subset)) colnames(design$X) else as.character(coordinate_subset)
  if (!length(coordinates) || anyNA(coordinates) || any(!nzchar(coordinates)) ||
      anyDuplicated(coordinates) || !all(coordinates %in% colnames(design$X))) {
    stop("WORKFLOW_COORDINATE_FAILED: all declared coordinates must remain estimable.")
  }
  if (!is.null(expected_coordinates)) {
    if (!setequal(coordinates, expected_coordinates)) {
      stop("WORKFLOW_COORDINATE_FAILED: deleting a subject changed the coordinate definitions.")
    }
    coordinates <- expected_coordinates
  }
  X <- design$X[, coordinates, drop = FALSE]
  if (isTRUE(random_slope)) {
    design$Z <- cbind(intercept = 1,
                      centered_time = metadata$cumulative_time - design$clinical$time_center)
  }
  n <- length(unique(metadata$subject_id))
  if (ncol(X) > floor(n / log(max(n, 3))) || qr(X)$rank != ncol(X)) {
    stop("LOWDIM_DIMENSION_FAILED: every workflow fit requires the declared dimension and full rank.")
  }
  list(
    y = log(metadata$c3), X = X, Z = design$Z, cluster_id = metadata$subject_id,
    preprocessing = list(
      clinical = design$clinical, module_rule = design$module_rule,
      gene_preprocessing = design$gene_preprocessing, module_audit = design$module_audit,
      sample_ids = rownames(metadata), subject_ids = unique(as.character(metadata$subject_id)),
      complete_design_coordinates = colnames(design$X), fitted_coordinates = coordinates,
      parameterisation = "sample_standardised_modules_and_sample_prepared_clinical_design",
      module_standardisation = "mean_and_sample_sd_recomputed_from_each_fitted_sample",
      coordinate_subset = coordinate_subset, random_slope = random_slope,
      coordinate_scope = if (is.null(coordinate_subset)) "complete_prespecified_design" else
        "declared_research_coordinate_subset"
    )
  )
}

lowdim_profile_workflow_diagnostics_v2 <- function(stage, omitted_subject, h) {
  fit <- stage$result$fit
  data.frame(
    omitted_subject = omitted_subject, n_clusters = length(stage$preprocessing$subject_ids),
    n_visits = length(stage$preprocessing$sample_ids),
    n_coordinates = length(stage$preprocessing$fitted_coordinates), h = h, h_half = h / 2,
    converged = fit$converged, refit_gradient = fit$post_refit_gradient_max,
    nuisance_gradient = fit$components$max_nuisance_gradient,
    profile_identity_error = fit$components$profile_identity_error,
    hessian_min_eigenvalue = fit$lowdim_diagnostics$minimum_eigenvalue,
    hessian_condition_number = fit$lowdim_diagnostics$condition_number,
    inverse_defect = fit$lowdim_diagnostics$inverse_defect,
    stringsAsFactors = FALSE
  )
}

full_subject_profile_pair_workflow_jackknife_v2 <- function(
    metadata, expression_probe, annotation, module_signatures, module_score_members,
    module_rule, tau, h, lambda_gamma = 1, include_sledai = FALSE,
    control = list(), coordinate_subset = NULL, checkpoint_dir = NULL,
    random_slope = FALSE, checkpoint_storage = "complete") {
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
  if (!checkpoint_storage %in% c("complete", "compact")) stop("Unknown workflow checkpoint storage.")
  required <- c("subject_id", "c3")
  if (!all(required %in% names(metadata)) || !nrow(metadata) ||
      anyNA(metadata$subject_id) || any(!nzchar(as.character(metadata$subject_id))) ||
      any(!is.finite(metadata$c3)) || any(metadata$c3 <= 0) ||
      is.null(rownames(metadata)) || anyDuplicated(rownames(metadata)) ||
      !identical(colnames(expression_probe), rownames(metadata))) {
    stop("Invalid observed GSE65391 workflow records or expression identifiers.")
  }
  if (length(h) != 1L || !is.finite(h) || h <= 0 ||
      length(tau) != 1L || !is.finite(tau) || tau <= 0 || tau >= 1) {
    stop("The complete workflow requires one valid fixed bandwidth and quantile.")
  }
  if ("log_c3" %in% names(metadata) &&
      (any(!is.finite(metadata$log_c3)) || max(abs(metadata$log_c3 - log(metadata$c3))) > 1e-12)) {
    stop("Recorded log C3 differs from the positive observed C3 values.")
  }
  expression_probe <- subset_expression_to_frozen_modules_v2(
    expression_probe, annotation, module_score_members
  )
  subjects <- unique(as.character(metadata$subject_id))
  n <- length(subjects)
  if (n < 4L) stop("At least four subjects are required for a complete workflow jackknife.")
  if (!is.null(checkpoint_dir)) {
    if (dir.exists(checkpoint_dir) || file.exists(checkpoint_dir)) {
      stop("WORKFLOW_OUTPUT_EXISTS: choose a new versioned checkpoint directory.")
    }
    if (!dir.create(checkpoint_dir, recursive = TRUE)) stop("The workflow checkpoint directory could not be created.")
  }
  fit_stage <- function(md, expr, omitted_subject, expected_coordinates, file_id) {
    prepared <- prepare_gse65391_profile_workflow_stage_v2(
      md, expr, annotation, module_signatures, module_score_members, module_rule,
      include_sledai, coordinate_subset, expected_coordinates, random_slope
    )
    if (!is.null(checkpoint_dir)) {
      saved_preprocessing <- if (checkpoint_storage == "compact") {
        compact_application_preprocessing_v2(prepared$preprocessing)
      } else prepared$preprocessing
      saveRDS(saved_preprocessing,
              file.path(checkpoint_dir, paste0(file_id, "_preprocessing.rds")))
    }
    result <- fit_lowdim_profile_pair_v2(
      prepared$y, prepared$X, prepared$Z, prepared$cluster_id, tau, h,
      lambda_gamma = lambda_gamma, control = control
    )
    stage <- list(preprocessing = prepared$preprocessing, result = result,
                  omitted_subject = omitted_subject)
    stage$diagnostics <- lowdim_profile_workflow_diagnostics_v2(stage, omitted_subject, h)
    if (!is.null(checkpoint_dir)) {
      saved_stage <- if (checkpoint_storage == "compact") {
        compact_application_workflow_stage_v2(stage, paste0(file_id, "_preprocessing.rds"))
      } else stage
      saveRDS(saved_stage, file.path(checkpoint_dir, paste0(file_id, ".rds")))
    }
    stage
  }
  full <- fit_stage(metadata, expression_probe, NA_character_, NULL, "full_sample")
  coordinates <- full$preprocessing$fitted_coordinates
  estimates <- matrix(NA_real_, n, length(coordinates), dimnames = list(subjects, coordinates))
  leave_one <- setNames(vector("list", n), subjects)
  for (i in seq_along(subjects)) {
    keep <- as.character(metadata$subject_id) != subjects[i]
    leave_one[[i]] <- fit_stage(
      metadata[keep, , drop = FALSE], expression_probe[, keep, drop = FALSE],
      subjects[i], coordinates, sprintf("delete_subject_%04d", i)
    )
    beta <- leave_one[[i]]$result$fit$beta_pair
    if (!identical(names(beta), coordinates)) {
      stop("WORKFLOW_COORDINATE_FAILED: paired estimates changed their coordinate names.")
    }
    estimates[i, ] <- beta
  }
  if (nrow(estimates) != n || any(!is.finite(estimates)) ||
      !identical(names(leave_one), subjects)) {
    stop("JACKKNIFE_INCOMPLETE: every subject requires a complete preprocessing and paired fit.")
  }
  centered <- sweep(estimates, 2, colMeans(estimates))
  diagnostics <- do.call(rbind, lapply(leave_one, function(stage) stage$diagnostics))
  if (!is.null(checkpoint_dir)) {
    utils::write.csv(diagnostics, file.path(checkpoint_dir, "leave_one_diagnostics.csv"), row.names = FALSE)
  }
  list(
    fit = full$result$fit, table = full$result$table, precision = full$result$precision,
    full_sample = full, leave_one = leave_one, estimates = estimates,
    se = sqrt((n - 1) / n * colSums(centered^2)), diagnostics = diagnostics,
    n_clusters = n, bandwidth_rule = "fixed_full_sample_bandwidth_in_leave_one_fits",
    jackknife_scope = "complete_sample_standardised_workflow",
    parameterisation = "sample_standardised_modules_and_sample_prepared_clinical_design",
    coordinate_scope = full$preprocessing$coordinate_scope,
    coordinate_subset = coordinate_subset, checkpoint_dir = checkpoint_dir,
    random_slope = random_slope, checkpoint_storage = checkpoint_storage
  )
}

lowdim_profile_pair_workflow_inference_v2 <- function(
    metadata, expression_probe, annotation, module_signatures, module_score_members,
    module_rule, tau, h, lambda_gamma = 1, include_sledai = FALSE,
    ci_level = 0.95, control = list(), coordinate_subset = NULL, checkpoint_dir = NULL,
    random_slope = FALSE, checkpoint_storage = "complete") {
  if (length(ci_level) != 1L || !is.finite(ci_level) || ci_level <= 0 || ci_level >= 1) {
    stop("Invalid confidence level.")
  }
  jack <- full_subject_profile_pair_workflow_jackknife_v2(
    metadata, expression_probe, annotation, module_signatures, module_score_members,
    module_rule, tau, h, lambda_gamma, include_sledai, control, coordinate_subset, checkpoint_dir,
    random_slope, checkpoint_storage
  )
  beta <- as.numeric(jack$fit$beta_pair)
  critical <- stats::qnorm((1 + ci_level) / 2)
  sandwich_se <- jack$table$se
  table <- data.frame(
    coordinate = colnames(jack$estimates), beta_hat = beta,
    beta_refit = as.numeric(jack$fit$beta_refit), score_source = "bandwidth_pair",
    score_definition = jack$fit$score_definition, method_version = jack$fit$method_version,
    target_type = "unsmoothed_regularised_profile", h = h, h_half = h / 2,
    lambda_gamma = lambda_gamma, tau = tau, ci_level = ci_level, n_clusters = jack$n_clusters,
    random_slope = random_slope,
    parameterisation = jack$parameterisation, coordinate_scope = jack$coordinate_scope,
    sandwich_scope = "profile_score_given_full_sample_preprocessing",
    sandwich_se = sandwich_se, sandwich_lower = beta - critical * sandwich_se,
    sandwich_upper = beta + critical * sandwich_se, jackknife_scope = jack$jackknife_scope,
    jackknife_se = jack$se, jackknife_lower = beta - critical * jack$se,
    jackknife_upper = beta + critical * jack$se,
    jackknife_completed = nrow(jack$estimates), stringsAsFactors = FALSE
  )
  list(table = table, fit = jack$fit, jackknife = jack, precision = jack$precision,
       preprocessing = jack$full_sample$preprocessing)
}
