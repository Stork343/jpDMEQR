validate_gse65391_subject_folds_v2 <- function(metadata, folds, k = 5L, repetitions = 5L) {
  subjects <- unique(as.character(metadata$subject_id))
  if (!all(c("subject_id", "repetition", "fold") %in% names(folds)) ||
      anyNA(folds[c("subject_id", "repetition", "fold")]) ||
      anyDuplicated(folds[c("subject_id", "repetition")]) ||
      !setequal(folds$repetition, seq_len(repetitions))) {
    stop("SUBJECT_FOLD_FAILED: incomplete or duplicate subject assignments.")
  }
  summary <- list()
  for (r in seq_len(repetitions)) {
    rows <- folds[folds$repetition == r, , drop = FALSE]
    if (!setequal(as.character(rows$subject_id), subjects) || nrow(rows) != length(subjects) ||
        !setequal(rows$fold, seq_len(k))) {
      stop("SUBJECT_FOLD_FAILED: every subject must appear once in every repetition.")
    }
    for (f in seq_len(k)) {
      test_subjects <- as.character(rows$subject_id[rows$fold == f])
      train_subjects <- setdiff(subjects, test_subjects)
      test <- as.character(metadata$subject_id) %in% test_subjects
      if (length(intersect(train_subjects, test_subjects)) || !any(test) || all(test)) {
        stop("SUBJECT_FOLD_FAILED: subject leakage or an empty role was found.")
      }
      summary[[length(summary) + 1L]] <- data.frame(
        repetition = r, fold = f, training_subjects = length(train_subjects),
        test_subjects = length(test_subjects), training_visits = sum(!test), test_visits = sum(test),
        overlap_subjects = 0L, stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, summary)
}

application_module_simultaneous_intervals_v2 <- function(table, family_size = 30L, alpha = 0.05) {
  table$coordinate_type <- ifelse(grepl("^module__", table$coordinate), "module", "clinical")
  table$simultaneous_family_size <- ifelse(table$coordinate_type == "module", family_size, NA_integer_)
  critical <- stats::qnorm(1 - alpha / (2 * family_size))
  module <- table$coordinate_type == "module"
  for (source in c("sandwich", "jackknife")) {
    lower <- upper <- rep(NA_real_, nrow(table))
    lower[module] <- table$beta_hat[module] - critical * table[[paste0(source, "_se")]][module]
    upper[module] <- table$beta_hat[module] + critical * table[[paste0(source, "_se")]][module]
    table[[paste0(source, "_simultaneous_lower")]] <- lower
    table[[paste0(source, "_simultaneous_upper")]] <- upper
  }
  table
}

prepare_gse65391_prediction_design_v2 <- function(
    metadata, expression_probe, annotation, module_signatures, module_score_members,
    module_rule, training_rows, include_sledai = FALSE, filter_size = 5000L,
    screening_dimension = NULL, screening_tau = NULL, random_slope = FALSE) {
  design <- prepare_gse65391_lowdim_design_v2(
    metadata, expression_probe, annotation, module_signatures, training_rows,
    module_rule, include_sledai, module_score_members
  )
  train <- if (is.logical(training_rows)) which(training_rows) else as.integer(training_rows)
  filtered <- remove_probe_rows_v2(expression_probe, annotation)
  collapsed <- collapse_probes_highest_training_mad_v2(
    filtered$expression, filtered$annotation$gene_symbol, rownames(metadata)[train]
  )
  selected <- filter_top_training_mad_v2(collapsed$expression, rownames(metadata)[train], filter_size)
  scaled <- scale_from_training_v2(t(selected$expression), train)
  genes <- scaled$X
  if (!is.null(screening_dimension)) {
    if (is.null(screening_tau) || !is.finite(screening_tau) ||
        screening_tau <= 0 || screening_tau >= 1) stop("Supply the fixed screening quantile.")
    cluster <- as.character(metadata$subject_id[train])
    residual_score <- screening_tau - as.numeric(metadata$log_c3[train] <=
                                                   stats::quantile(metadata$log_c3[train], screening_tau,
                                                                   names = FALSE, type = 7))
    contributions <- genes[train, , drop = FALSE] * residual_score
    by_subject <- rowsum(contributions, cluster, reorder = FALSE)
    by_subject <- by_subject / as.numeric(table(factor(cluster, levels = rownames(by_subject))))
    score <- abs(colMeans(by_subject))
    selected_genes <- order(score, decreasing = TRUE)[seq_len(min(as.integer(screening_dimension), ncol(genes)))]
    genes <- genes[, selected_genes, drop = FALSE]
    screening <- data.frame(gene_symbol = colnames(scaled$X), score = score,
                             selected = seq_len(ncol(scaled$X)) %in% selected_genes,
                             stringsAsFactors = FALSE)
  } else screening <- NULL
  colnames(genes) <- paste0("gene__", colnames(genes))
  design$X <- cbind(design$X, genes)
  design$penalty_factor <- c(rep(0, ncol(design$X) - ncol(genes)), rep(1, ncol(genes)))
  design$Z <- if (isTRUE(random_slope)) {
    cbind(intercept = 1, centered_time = metadata$cumulative_time - design$clinical$time_center)
  } else design$Z
  design$highdim_preprocessing <- list(
    filter_size = filter_size, screening_dimension = screening_dimension,
    screening_tau = screening_tau, screening = screening,
    gene_ids = colnames(scaled$X), center = scaled$center, scale = scaled$scale,
    selected_genes = colnames(genes),
    gene_probe_map = data.frame(gene_symbol = rownames(collapsed$expression),
                                probe_id = rownames(filtered$expression)[collapsed$selected_probe_rows],
                                stringsAsFactors = FALSE),
    training_subjects = unique(as.character(metadata$subject_id[train]))
  )
  design
}
