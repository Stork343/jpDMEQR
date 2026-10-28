compact_application_preprocessing_v2 <- function(preprocessing) {
  result <- preprocessing
  result$clinical$X <- NULL
  environment(result$clinical$formula) <- baseenv()
  result
}

compact_application_component_diagnostics_v2 <- function(component) {
  fields <- c("objective", "score", "profile_identity_error", "identity_checked",
              "cluster_identity_checked", "cluster_identity_errors", "matrix_backend",
              "nuisance_diagnostics", "nuisance_converged", "max_nuisance_gradient",
              "nuisance_iterations", "Lambda", "n_clusters", "total_rows", "p", "q",
              "tau", "h", "hessian_min_eigenvalue", "bandwidths", "bandwidth_weights",
              "score_definition")
  result <- component[intersect(fields, names(component))]
  if (!is.null(component$components_h)) {
    result$components_h <- compact_application_component_diagnostics_v2(component$components_h)
    result$components_half_h <- compact_application_component_diagnostics_v2(component$components_half_h)
  }
  result
}

compact_application_workflow_stage_v2 <- function(stage, preprocessing_file) {
  fit <- stage$result$fit
  fit$components <- compact_application_component_diagnostics_v2(fit$components)
  fit$components_first_stage <- NULL
  fit$call <- NULL
  list(result = list(fit = fit), diagnostics = stage$diagnostics,
       omitted_subject = stage$omitted_subject, preprocessing_file = preprocessing_file,
       storage_mode = "complete_diagnostics_and_preprocessing_reference")
}

save_compact_application_workflow_v2 <- function(result, path, extra = list()) {
  if (file.exists(path)) stop("The compact workflow artifact already exists.")
  directory <- result$jackknife$checkpoint_dir
  if (is.null(directory) || result$jackknife$checkpoint_storage != "compact") {
    stop("Compact workflow storage requires compact checkpoints from every fitted sample.")
  }
  subjects <- rownames(result$jackknife$estimates)
  index <- data.frame(
    omitted_subject = subjects,
    fit_file = sprintf("delete_subject_%04d.rds", seq_along(subjects)),
    preprocessing_file = sprintf("delete_subject_%04d_preprocessing.rds", seq_along(subjects)),
    stringsAsFactors = FALSE
  )
  required <- c("full_sample.rds", "full_sample_preprocessing.rds", index$fit_file, index$preprocessing_file)
  if (any(!file.exists(file.path(directory, required)))) stop("A complete workflow checkpoint is missing.")
  jack_fields <- setdiff(names(result$jackknife), c("fit", "table", "precision", "full_sample", "leave_one", "checkpoint_dir"))
  artifact <- list(
    storage_mode = "complete_workflow_compact_v1", fit = result$fit,
    table = result$table, precision = result$precision,
    jackknife = result$jackknife[jack_fields], leave_one_index = index,
    checkpoint_directory = basename(directory), extra = extra
  )
  saveRDS(artifact, path)
  invisible(path)
}

read_compact_application_workflow_v2 <- function(path) {
  artifact <- readRDS(path)
  if (!identical(artifact$storage_mode, "complete_workflow_compact_v1")) {
    stop("The artifact is not a complete compact workflow.")
  }
  directory <- file.path(dirname(path), artifact$checkpoint_directory)
  read_stage <- function(fit_file) {
    stage <- readRDS(file.path(directory, fit_file))
    stage$preprocessing <- readRDS(file.path(directory, stage$preprocessing_file))
    stage
  }
  index <- artifact$leave_one_index
  leave <- lapply(index$fit_file, read_stage)
  names(leave) <- index$omitted_subject
  estimates <- do.call(rbind, lapply(leave, function(stage) stage$result$fit$beta_pair))
  rownames(estimates) <- index$omitted_subject
  if (!identical(estimates, artifact$jackknife$estimates) ||
      any(!vapply(leave, function(stage) isTRUE(stage$result$fit$converged), logical(1)))) {
    stop("Compact leave-one coefficients or convergence records changed.")
  }
  n <- nrow(estimates)
  centered <- sweep(estimates, 2, colMeans(estimates))
  jack <- artifact$jackknife
  jack$se <- sqrt((n - 1) / n * colSums(centered^2))
  jack$leave_one <- leave
  jack$checkpoint_dir <- directory
  jack$fit <- artifact$fit
  jack$precision <- artifact$precision
  jack$table <- profile_bandwidth_pair_wald_v2(
    artifact$fit$components, artifact$fit$beta_refit, artifact$precision, seq_len(ncol(estimates))
  )
  jack$full_sample <- read_stage("full_sample.rds")
  jack$full_sample$result <- list(fit = artifact$fit, table = jack$table, precision = artifact$precision)
  table <- artifact$table
  critical <- stats::qnorm((1 + table$ci_level[1]) / 2)
  table$beta_hat <- as.numeric(jack$table$beta_tilde)
  table$beta_refit <- as.numeric(artifact$fit$beta_refit)
  table$sandwich_se <- jack$table$se
  table$sandwich_lower <- table$beta_hat - critical * table$sandwich_se
  table$sandwich_upper <- table$beta_hat + critical * table$sandwich_se
  table$jackknife_se <- jack$se
  table$jackknife_lower <- table$beta_hat - critical * table$jackknife_se
  table$jackknife_upper <- table$beta_hat + critical * table$jackknife_se
  table$jackknife_completed <- n
  if (!isTRUE(all.equal(table, artifact$table, tolerance = 0, check.attributes = TRUE))) {
    stop("Compact readback did not reproduce the complete inference table.")
  }
  list(table = table, fit = artifact$fit, jackknife = jack,
       precision = artifact$precision, preprocessing = jack$full_sample$preprocessing,
       storage_extra = artifact$extra)
}

application_matched_coordinate_metrics_v2 <- function(result, population, target, case, repeat_id) {
  coordinates <- names(population$beta0)
  n <- result$jackknife$n_clusters
  if (!identical(names(result$fit$beta_pair), coordinates)) stop("The matched categorical coordinate design changed.")
  full_map <- application_population_unit_map_v2(result$preprocessing, population$preprocessing, coordinates)
  beta <- as.numeric(full_map$matrix %*% result$fit$beta_pair)
  G <- result$fit$components$cluster_scores
  influence <- sweep(G, 2, colMeans(G)) %*% result$precision
  V <- full_map$matrix %*% (crossprod(influence) / n^2) %*% t(full_map$matrix)
  sandwich_se <- sqrt(diag(V))
  estimates <- t(vapply(result$jackknife$leave_one, function(stage) {
    if (!identical(names(stage$result$fit$beta_pair), coordinates)) stop("A matched leave-one coordinate design changed.")
    map <- application_population_unit_map_v2(stage$preprocessing, population$preprocessing, coordinates)
    as.numeric(map$matrix %*% stage$result$fit$beta_pair)
  }, numeric(length(coordinates))))
  colnames(estimates) <- coordinates
  centered <- sweep(estimates, 2, colMeans(estimates))
  jackknife_se <- sqrt((n - 1) / n * colSums(centered^2))
  target <- target[match(coordinates, target$coordinate), ]
  critical <- stats::qnorm(.975)
  table <- data.frame(
    case_id = case$case_id, model = case$model, q = case$q, error_family = case$error_family,
    tau = case$tau, replicate = repeat_id, coordinate = coordinates, beta_hat = beta,
    target = target$target, target_mcse = target$target_mcse, target_source = target$target_source,
    sandwich_se = sandwich_se, jackknife_se = jackknife_se,
    sandwich_lower = beta - critical * sandwich_se, sandwich_upper = beta + critical * sandwich_se,
    jackknife_lower = beta - critical * jackknife_se, jackknife_upper = beta + critical * jackknife_se,
    changed_gene_probes = full_map$changed_gene_probes,
    exact_affine_measurement_map = full_map$exact_affine_measurement_map,
    reporting_units = "fixed_population_preprocessing_units", status = "ok", stringsAsFactors = FALSE
  )
  table$sandwich_covered <- table$sandwich_lower <= table$target & table$sandwich_upper >= table$target
  table$jackknife_covered <- table$jackknife_lower <= table$target & table$jackknife_upper >= table$target
  table <- application_module_simultaneous_intervals_v2(table)
  table$sandwich_simultaneous_covered <- table$sandwich_simultaneous_lower <= table$target &
    table$sandwich_simultaneous_upper >= table$target
  table$jackknife_simultaneous_covered <- table$jackknife_simultaneous_lower <= table$target &
    table$jackknife_simultaneous_upper >= table$target
  list(table = table, population_unit_map = full_map, mapped_jackknife = estimates)
}
