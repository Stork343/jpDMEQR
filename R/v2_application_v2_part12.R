validate_exploratory_subject_roles_v2 <- function(metadata, roles, splits = 10L) {
  required <- c("split", "swapped", "subject_id", "role")
  if (!all(required %in% names(roles)) || anyNA(roles[required]) ||
      !setequal(unique(roles$split), seq_len(splits)) ||
      !setequal(unique(roles$swapped), c(FALSE, TRUE)) ||
      !all(roles$role %in% c("selection", "inference")) ||
      anyDuplicated(roles[c("split", "swapped", "subject_id")])) {
    stop("The complete frozen exploratory subject-role table is required.")
  }
  subjects <- unique(as.character(metadata$subject_id))
  output <- list()
  for (split in seq_len(splits)) {
    original <- roles[roles$split == split & !roles$swapped, , drop = FALSE]
    reversed <- roles[roles$split == split & roles$swapped, , drop = FALSE]
    if (!setequal(original$subject_id, subjects) || !setequal(reversed$subject_id, subjects) ||
        nrow(original) != length(subjects) || nrow(reversed) != length(subjects)) {
      stop("Every exploratory role assignment must contain every eligible subject once.")
    }
    if (any(original$role == reversed$role[match(original$subject_id, reversed$subject_id)])) {
      stop("The swapped exploratory assignment must exchange every subject's role.")
    }
    for (swapped in c(FALSE, TRUE)) {
      assignment <- roles[roles$split == split & roles$swapped == swapped, , drop = FALSE]
      selection <- as.character(assignment$subject_id[assignment$role == "selection"])
      inference <- as.character(assignment$subject_id[assignment$role == "inference"])
      if (!length(selection) || !length(inference) ||
          abs(length(selection) - length(inference)) > 1L || length(intersect(selection, inference))) {
        stop("Exploratory roles require disjoint subject halves.")
      }
      output[[length(output) + 1L]] <- data.frame(
        split = split, swapped = swapped,
        role_id = sprintf("split%02d_swap%d", split, as.integer(swapped)),
        selection_subjects = length(selection), inference_subjects = length(inference),
        selection_visits = sum(metadata$subject_id %in% selection),
        inference_visits = sum(metadata$subject_id %in% inference), stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, output)
}

prepare_exploratory_gene_split_v2 <- function(metadata, expression_probe, annotation,
                                               assignment, filter_size = 5000L) {
  if (!all(c("subject_id", "role") %in% names(assignment)) ||
      anyDuplicated(assignment$subject_id) ||
      !setequal(assignment$subject_id, unique(metadata$subject_id)) ||
      !identical(colnames(expression_probe), rownames(metadata))) {
    stop("The exploratory expression, metadata and frozen assignment must agree.")
  }
  train <- as.character(metadata$subject_id) %in%
    as.character(assignment$subject_id[assignment$role == "selection"])
  test <- as.character(metadata$subject_id) %in%
    as.character(assignment$subject_id[assignment$role == "inference"])
  if (any(train & test) || any(!(train | test)) || !any(train) || !any(test)) {
    stop("The selection and inference subjects must partition the eligible cohort.")
  }
  filtered <- remove_probe_rows_v2(expression_probe, annotation)
  collapsed <- collapse_probes_highest_training_mad_v2(
    filtered$expression, filtered$annotation$gene_symbol, rownames(metadata)[train]
  )
  selected <- filter_top_training_mad_v2(collapsed$expression, rownames(metadata)[train], filter_size)
  scaled <- scale_from_training_v2(t(selected$expression), train)
  if (ncol(scaled$X) < 3L || any(!is.finite(scaled$X))) {
    stop("At least three finite training-standardised genes are required.")
  }
  clinical <- lapply(c("A", "B"), function(model) {
    prepare_clinical_design_split_v2(metadata, train, include_sledai = model == "B")
  })
  names(clinical) <- c("A", "B")
  map <- data.frame(
    gene_id = rownames(collapsed$expression),
    probe_id = rownames(filtered$expression)[collapsed$selected_probe_rows],
    stringsAsFactors = FALSE
  )
  list(genes = scaled$X, clinical = clinical, selection_rows = which(train),
       inference_rows = which(test), selection_subjects = unique(as.character(metadata$subject_id[train])),
       inference_subjects = unique(as.character(metadata$subject_id[test])),
       preprocessing = list(filter_size = as.integer(filter_size), gene_ids = colnames(scaled$X),
                            center = scaled$center, scale = scaled$scale,
                            probe_map = map, training_samples = rownames(metadata)[train]),
       assignment = assignment, sample_ids = rownames(metadata),
       preprocessing_scope = "selection_subjects_only")
}

exploratory_gene_ranking_v2 <- function(beta_l1, genes, y, cluster_id, tau,
                                        support_threshold = 1e-10) {
  genes <- assert_numeric_matrix_v2(genes, "genes")
  ids <- colnames(genes)
  coordinates <- paste0("gene__", ids)
  if (is.null(ids) || anyDuplicated(ids) || anyNA(match(coordinates, names(beta_l1))) ||
      length(y) != nrow(genes) || length(cluster_id) != nrow(genes) ||
      any(!is.finite(y)) || any(!is.finite(beta_l1[coordinates]))) {
    stop("The complete fitted gene coefficients and selection-role observations are required.")
  }
  cutoff <- stats::quantile(y, tau, names = FALSE, type = 7)
  contributions <- genes * (tau - as.numeric(y <= cutoff))
  by_subject <- rowsum(contributions, as.character(cluster_id), reorder = FALSE)
  counts <- as.numeric(table(factor(as.character(cluster_id), levels = rownames(by_subject))))
  score <- colMeans(by_subject / counts)
  table <- data.frame(gene_id = ids, coordinate = coordinates,
                      beta_l1 = as.numeric(beta_l1[coordinates]),
                      quantile_score = as.numeric(score),
                      L1_nonzero = abs(beta_l1[coordinates]) > support_threshold,
                      stringsAsFactors = FALSE)
  order_id <- order(-abs(table$beta_l1), -abs(table$quantile_score), table$gene_id,
                    method = "radix")
  table <- table[order_id, , drop = FALSE]
  rownames(table) <- NULL
  table$selection_rank <- seq_len(nrow(table))
  table$selected <- table$selection_rank <= 3L
  table$support_threshold <- support_threshold
  table
}

fit_exploratory_gene_role_v2 <- function(metadata, prepared, model, tau, checkpoint_dir,
                                          source_inputs, initial_control = list(),
                                          inference_control = list()) {
  if (!model %in% c("A", "B") || !tau %in% c(.25, .5, .75) ||
      !identical(rownames(metadata), prepared$sample_ids) || dir.exists(checkpoint_dir)) {
    stop("A declared model, quantile, matching preparation and new checkpoint directory are required.")
  }
  dir.create(checkpoint_dir, recursive = TRUE)
  selection_rows <- prepared$selection_rows
  inference_rows <- prepared$inference_rows
  clinical <- prepared$clinical[[model]]$X
  genes <- prepared$genes
  colnames(genes) <- paste0("gene__", colnames(genes))
  X <- cbind(clinical, genes)
  y <- as.numeric(metadata$log_c3)
  if (any(!is.finite(y)) || any(!is.finite(metadata$c3)) || any(metadata$c3 <= 0) ||
      max(abs(y - log(metadata$c3))) > 1e-12) stop("The frozen positive C3 outcome is required.")
  Z <- matrix(1, nrow(metadata), 1L, dimnames = list(rownames(metadata), "intercept"))
  n_inference <- length(unique(metadata$subject_id[inference_rows]))
  request <- list(model = model, tau = tau, source_inputs = source_inputs,
                  selection_subjects = prepared$selection_subjects,
                  inference_subjects = prepared$inference_subjects,
                  initial_control = initial_control, inference_control = inference_control,
                  h_inf = n_inference^(-3 / 10), lambda_gamma = 1,
                  mean_score_multiplier = 1, lambda_multiplier = 1,
                  ranking_rule = "descending_abs_beta_l1_then_abs_equal_subject_quantile_score_then_gene_id",
                  target_scope = "exploratory_role_specific_selected_design")
  saveRDS(request, file.path(checkpoint_dir, "request.rds"))
  stage <- "selection_fit"
  ranking <- NULL
  withCallingHandlers({
    initial <- fit_application_profile_initial_v2(
      y[selection_rows], X[selection_rows, , drop = FALSE], Z[selection_rows, , drop = FALSE],
      metadata$subject_id[selection_rows], tau,
      penalty_factor = c(rep(0, ncol(clinical)), rep(1, ncol(genes))), control = initial_control
    )
    ranking <- exploratory_gene_ranking_v2(
      initial$beta_l1, prepared$genes[selection_rows, , drop = FALSE], y[selection_rows],
      metadata$subject_id[selection_rows], tau
    )
    selected <- ranking$coordinate[ranking$selected]
    if (length(selected) != 3L) stop("The declared exploratory design requires exactly three genes.")
    reduced <- cbind(clinical, genes[, selected, drop = FALSE])
    if (ncol(reduced) > 15L || ncol(reduced) > floor(n_inference / log(n_inference)) ||
        qr(reduced[inference_rows, , drop = FALSE])$rank != ncol(reduced)) {
      stop("The independent inference design failed its dimension or rank requirement.")
    }
    selected_design_file <- file.path(checkpoint_dir, "selected_design.rds")
    saveRDS(list(X_selection = reduced[selection_rows, , drop = FALSE],
                 X_inference = reduced[inference_rows, , drop = FALSE],
                 y_selection = y[selection_rows], y_inference = y[inference_rows],
                 Z_inference = Z[inference_rows, , drop = FALSE],
                 cluster_inference = metadata$subject_id[inference_rows],
                 selected = ranking[ranking$selected, , drop = FALSE],
                 clinical = prepared$clinical[[model]], preprocessing = prepared$preprocessing,
                 source_inputs = source_inputs), selected_design_file)
    saveRDS(ranking, file.path(checkpoint_dir, "gene_ranking.rds"))
    initial$components <- NULL
    saveRDS(initial, file.path(checkpoint_dir, "selection_fit.rds"))
    rm(X, initial, genes)
    stage <- "independent_inference"
    inference <- fit_lowdim_profile_pair_v2(
      y[inference_rows], reduced[inference_rows, , drop = FALSE], Z[inference_rows, , drop = FALSE],
      metadata$subject_id[inference_rows], tau, h = request$h_inf,
      lambda_gamma = 1, control = inference_control
    )
    saveRDS(inference, file.path(checkpoint_dir, "inference_fit.rds"))
    table <- inference$table
    table$coordinate <- colnames(reduced)
    table <- table[table$coordinate %in% selected, , drop = FALSE]
    table$gene_id <- sub("^gene__", "", table$coordinate)
    table$ci_level <- .95
    table$status <- "ok"
    table$failure_stage <- NA_character_
    table$n_clusters <- n_inference
    table$selected_design_file <- normalizePath(selected_design_file)
    table$target_scope <- request$target_scope
    table$method_version <- inference$fit$method_version
    table$score_definition <- inference$fit$score_definition
    utils::write.csv(table, file.path(checkpoint_dir, "intervals.csv"), row.names = FALSE)
    result <- list(status = "ok", stage = "completed", ranking = ranking, table = table,
                   n_clusters = n_inference, selected_design_file = normalizePath(selected_design_file),
                   lowdim_diagnostics = inference$fit$lowdim_diagnostics,
                   refit_gradient = inference$fit$post_refit_gradient_max,
                   nuisance_gradient = inference$fit$components$max_nuisance_gradient,
                   profile_identity_error = inference$fit$components$profile_identity_error)
    saveRDS(result, file.path(checkpoint_dir, "result.rds"))
    result
  }, error = function(condition) {
    result <- list(status = "failed", stage = condition$failure_stage %||% stage,
                   message = conditionMessage(condition), condition = condition,
                   ranking = ranking, n_clusters = n_inference,
                   selected_design_file = if (file.exists(file.path(checkpoint_dir, "selected_design.rds")))
                     normalizePath(file.path(checkpoint_dir, "selected_design.rds")) else NA_character_)
    saveRDS(result, file.path(checkpoint_dir, "result.rds"))
  })
}

summarise_exploratory_gene_roles_v2 <- function(jobs, results) {
  if (nrow(jobs) != length(results) || anyDuplicated(jobs$job_id)) stop("Exploratory jobs and results disagree.")
  selection <- intervals <- metrics <- list()
  for (i in seq_len(nrow(jobs))) {
    job <- jobs[i, , drop = FALSE]
    result <- results[[i]]
    status <- if (is.null(result)) "not_attempted" else result$status
    metrics[[i]] <- cbind(job, status = status,
                          failure_stage = if (status == "failed") result$stage else NA_character_,
                          failure_message = if (status == "failed") result$message else NA_character_,
                          selection_completed = !is.null(result$ranking),
                          n_clusters = if (is.null(result)) NA_integer_ else result$n_clusters)
    if (is.null(result$ranking)) next
    selection[[length(selection) + 1L]] <- data.frame(
      job[rep(1L, nrow(result$ranking)), , drop = FALSE], result$ranking, row.names = NULL
    )
    selected <- result$ranking[result$ranking$selected, , drop = FALSE]
    table <- if (status == "ok") result$table else data.frame(
      coordinate = selected$coordinate, gene_id = selected$gene_id,
      beta_tilde = NA_real_, se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
      ci_level = .95, status = "failed", failure_stage = result$stage,
      n_clusters = result$n_clusters, selected_design_file = result$selected_design_file,
      target_scope = "exploratory_role_specific_selected_design", stringsAsFactors = FALSE
    )
    fields <- c("gene_id", "coordinate", "beta_tilde", "se", "ci_lower", "ci_upper",
                "ci_level", "status", "failure_stage", "n_clusters", "selected_design_file", "target_scope")
    intervals[[length(intervals) + 1L]] <- data.frame(
      job[rep(1L, nrow(table)), , drop = FALSE], table[, fields, drop = FALSE], row.names = NULL
    )
  }
  metrics <- do.call(rbind, metrics)
  selection <- if (length(selection)) do.call(rbind, selection) else data.frame()
  intervals <- if (length(intervals)) do.call(rbind, intervals) else data.frame()
  frequency <- coefficients <- list()
  if (nrow(selection)) {
    groups <- unique(selection[c("model", "tau")])
    for (i in seq_len(nrow(groups))) {
      group <- groups[i, , drop = FALSE]
      ranks <- selection[selection$model == group$model & selection$tau == group$tau, , drop = FALSE]
      attempts <- metrics[metrics$model == group$model & metrics$tau == group$tau, , drop = FALSE]
      scheduled <- nrow(attempts)
      completed <- sum(attempts$selection_completed)
      for (gene in sort(unique(ranks$gene_id))) {
        rows <- ranks[ranks$gene_id == gene, , drop = FALSE]
        count <- sum(rows$selected)
        frequency[[length(frequency) + 1L]] <- data.frame(
          model = group$model, tau = group$tau, gene_id = gene, scheduled_roles = scheduled,
          selection_completed_roles = completed, eligible_roles = nrow(rows), selected_roles = count,
          L1_nonzero_selected_roles = sum(rows$selected & rows$L1_nonzero),
          selection_frequency = count / completed, scheduled_selection_fraction = count / scheduled,
          stringsAsFactors = FALSE
        )
        if (!count) next
        values <- intervals[intervals$model == group$model & intervals$tau == group$tau &
                              intervals$gene_id == gene & intervals$status == "ok", "beta_tilde"]
        quantiles <- if (length(values)) stats::quantile(values, c(.25, .5, .75), names = FALSE) else
          rep(NA_real_, 3L)
        coefficients[[length(coefficients) + 1L]] <- data.frame(
          model = group$model, tau = group$tau, gene_id = gene, n_selected_roles = count,
          n_inference_ok_roles = length(values), beta_tilde_median = quantiles[2],
          beta_tilde_q25 = quantiles[1], beta_tilde_q75 = quantiles[3],
          target_scope = "descriptive_selected_designs", stringsAsFactors = FALSE
        )
      }
    }
  }
  list(role_metrics = metrics, selection = selection, intervals = intervals,
       selection_frequency = if (length(frequency)) do.call(rbind, frequency) else data.frame(),
       coefficient_median_iqr = if (length(coefficients)) do.call(rbind, coefficients) else data.frame())
}
