audit_primary_outcome_v2 <- function(metadata,
                                     outcome = "c3",
                                     minimum_subjects = 100L,
                                     minimum_visits = 500L,
                                     minimum_repeated_subjects = 75L,
                                     maximum_missing_fraction = 0.30,
                                     minimum_distinct_values = 20L,
                                     maximum_modal_fraction = 0.20,
                                     minimum_visits_per_subject = 2L) {
  if (!outcome %in% names(metadata)) stop("Outcome field not found: ", outcome)
  sle <- metadata$group_type == "SLE" & !metadata$technical_replicate_flag
  y_all <- suppressWarnings(as.numeric(metadata[[outcome]]))
  eligible_visit <- sle & is.finite(y_all) & y_all > 0 & !is.na(metadata$subject_id)
  subjects_all <- unique(metadata$subject_id[sle & !is.na(metadata$subject_id)])
  subjects_eligible <- unique(metadata$subject_id[eligible_visit])
  visit_counts <- table(metadata$subject_id[eligible_visit])
  repeated_subjects <- names(visit_counts)[visit_counts >= minimum_visits_per_subject]
  final_visit <- eligible_visit & metadata$subject_id %in% repeated_subjects
  y <- y_all[final_visit]

  missing_fraction <- safe_fraction_v2(sum(sle & !is.finite(y_all)), sum(sle))
  modal_fraction <- if (length(y)) max(table(y)) / length(y) else NA_real_
  gates <- data.frame(
    gate = c(
      "minimum_subjects", "minimum_visits", "minimum_repeated_subjects",
      "maximum_missing_fraction", "positive_values", "minimum_distinct_values",
      "maximum_modal_fraction"
    ),
    observed = c(
      length(unique(metadata$subject_id[final_visit])),
      sum(final_visit),
      length(repeated_subjects),
      missing_fraction,
      if (length(y)) min(y) else NA_real_,
      length(unique(y)),
      modal_fraction
    ),
    threshold = c(
      minimum_subjects, minimum_visits, minimum_repeated_subjects,
      maximum_missing_fraction, 0, minimum_distinct_values, maximum_modal_fraction
    ),
    direction = c(">=", ">=", ">=", "<=", ">", ">=", "<="),
    pass = c(
      length(unique(metadata$subject_id[final_visit])) >= minimum_subjects,
      sum(final_visit) >= minimum_visits,
      length(repeated_subjects) >= minimum_repeated_subjects,
      is.finite(missing_fraction) && missing_fraction <= maximum_missing_fraction,
      length(y) > 0 && all(y > 0),
      length(unique(y)) >= minimum_distinct_values,
      is.finite(modal_fraction) && modal_fraction <= maximum_modal_fraction
    ),
    stringsAsFactors = FALSE
  )

  summary <- data.frame(
    source_sle_subjects = length(subjects_all),
    source_sle_visits = sum(sle),
    outcome_observed_subjects = length(subjects_eligible),
    outcome_observed_visits = sum(eligible_visit),
    repeated_eligible_subjects = length(repeated_subjects),
    final_subjects = length(unique(metadata$subject_id[final_visit])),
    final_visits = sum(final_visit),
    missing_fraction = missing_fraction,
    distinct_values = length(unique(y)),
    modal_fraction = modal_fraction,
    outcome_min = if (length(y)) min(y) else NA_real_,
    outcome_median = if (length(y)) stats::median(y) else NA_real_,
    outcome_max = if (length(y)) max(y) else NA_real_,
    all_gates_pass = all(gates$pass),
    stringsAsFactors = FALSE
  )

  list(
    gates = gates,
    summary = summary,
    eligible_visit = final_visit,
    eligible_subjects = repeated_subjects,
    outcome = y_all,
    transformed_outcome = ifelse(final_visit, log(y_all), NA_real_)
  )
}

cohort_flow_gse65391_v2 <- function(metadata, eligible_visit = NULL) {
  eligible_visit <- eligible_visit %||% rep(FALSE, nrow(metadata))
  steps <- data.frame(
    step = c(
      "GEO samples",
      "SLE biological samples",
      "Observed positive primary outcome",
      "Subjects with >=2 eligible visits",
      "Final eligible visits"
    ),
    samples = c(
      nrow(metadata),
      sum(metadata$group_type == "SLE" & !metadata$technical_replicate_flag),
      sum(metadata$group_type == "SLE" & is.finite(metadata$c3) & metadata$c3 > 0,
          na.rm = TRUE),
      sum(metadata$subject_id %in% unique(metadata$subject_id[eligible_visit])),
      sum(eligible_visit)
    ),
    subjects = c(
      length(unique(metadata$subject_id[!is.na(metadata$subject_id)])),
      length(unique(metadata$subject_id[metadata$group_type == "SLE" &
                                          !metadata$technical_replicate_flag])),
      length(unique(metadata$subject_id[metadata$group_type == "SLE" &
                                          is.finite(metadata$c3) & metadata$c3 > 0])),
      length(unique(metadata$subject_id[eligible_visit])),
      length(unique(metadata$subject_id[eligible_visit]))
    ),
    stringsAsFactors = FALSE
  )
  steps
}

make_repeated_subject_folds_v2 <- function(metadata,
                                           eligible_visit,
                                           k = 5L,
                                           repetitions = 5L,
                                           seed = 20260817L) {
  md <- metadata[eligible_visit, , drop = FALSE]
  subjects <- unique(md$subject_id)
  if (length(subjects) < k) stop("Too few subjects for requested folds.")
  batch_by_subject <- tapply(md$batch, md$subject_id, function(x) {
    x <- x[is.finite(x)]
    if (length(x)) as.character(names(sort(table(x), decreasing = TRUE))[1]) else "missing"
  })
  baseline_c3 <- tapply(md$c3, md$subject_id, function(x) {
    pos <- which(is.finite(x))
    if (length(pos)) x[pos[1]] else NA_real_
  })
  finite_baseline <- is.finite(baseline_c3)
  c3_group <- rep("missing", length(baseline_c3))
  names(c3_group) <- names(baseline_c3)
  if (sum(finite_baseline) >= 3L) {
    ranked <- rank(baseline_c3[finite_baseline], ties.method = "average")
    # Rank-based thirds avoid duplicated empirical-quantile breaks when C3 is discrete.
    c3_group[finite_baseline] <- as.character(
      pmin(3L, pmax(1L, ceiling(3 * ranked / length(ranked))))
    )
  }
  strata <- paste(batch_by_subject[subjects], c3_group[subjects], sep = "::")

  out <- vector("list", repetitions)
  for (r in seq_len(repetitions)) {
    set.seed(seed + r * 1009L)
    fold <- integer(length(subjects))
    for (st in unique(strata)) {
      idx <- which(strata == st)
      fold[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
    }
    out[[r]] <- data.frame(
      subject_id = subjects,
      repetition = r,
      fold = fold,
      stratum = strata,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

remove_probe_rows_v2 <- function(expression,
                                 annotation,
                                 probe_col = "probe_id",
                                 gene_col = "gene_symbol",
                                 control_col = "is_control",
                                 multimapping_col = "is_multimapping") {
  expression <- as.matrix(expression)
  annotation <- as.data.frame(annotation, stringsAsFactors = FALSE)
  if (!probe_col %in% names(annotation)) stop("Annotation probe column missing.")
  idx <- match(rownames(expression), annotation[[probe_col]])
  ann <- annotation[idx, , drop = FALSE]
  gene <- standardize_geo_missing_v2(ann[[gene_col]])
  keep <- !is.na(gene)
  if (control_col %in% names(ann)) keep <- keep & !as.logical(ann[[control_col]])
  if (multimapping_col %in% names(ann)) keep <- keep & !as.logical(ann[[multimapping_col]])
  list(expression = expression[keep, , drop = FALSE], annotation = ann[keep, , drop = FALSE],
       keep = keep)
}

collapse_probes_highest_training_mad_v2 <- function(expression,
                                                    gene_symbol,
                                                    training_samples) {
  expression <- as.matrix(expression)
  gene_symbol <- as.character(gene_symbol)
  train_idx <- match(training_samples, colnames(expression))
  if (anyNA(train_idx)) stop("Unknown training sample identifiers.")
  mad_probe <- apply(expression[, train_idx, drop = FALSE], 1, stats::mad,
                     constant = 1, na.rm = TRUE)
  by_gene <- split(seq_along(gene_symbol), gene_symbol)
  selected <- vapply(by_gene, function(idx) idx[which.max(mad_probe[idx])], integer(1))
  out <- expression[selected, , drop = FALSE]
  rownames(out) <- names(selected)
  list(expression = out, selected_probe_rows = selected, training_mad = mad_probe[selected])
}
