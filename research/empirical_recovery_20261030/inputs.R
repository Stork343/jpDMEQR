empirical_atomic_json <- function(value, path) {
  temporary <- paste0(path, ".", Sys.getpid(), ".writing")
  jsonlite::write_json(value, temporary, auto_unbox = TRUE, pretty = TRUE)
  if (!file.rename(temporary, path)) stop("Unable to publish the empirical checks.")
}

empirical_inputs <- function(load_expression = TRUE) {
  cohort_path <- "data/derived/GSE65391/audit_20261022_bootstrap_v1/analysis_cohort.rds"
  cohort <- readRDS(cohort_path)
  metadata <- cohort$eligible_metadata
  if (nrow(metadata) != 867L || length(unique(metadata$subject_id)) != 129L ||
      anyDuplicated(rownames(metadata)) || any(!is.finite(metadata$log_c3)) ||
      any(metadata$c3 <= 0) || max(abs(metadata$log_c3 - log(metadata$c3))) > 1e-12)
    stop("The frozen C3 cohort does not match its declared observations.")
  modules <- validate_gse65391_fit_inputs_v2(cohort$audit, "config/module_signatures.csv")
  members <- read_frozen_module_score_members_v2("config/module_score_members_20261022.csv")
  rule <- yaml::read_yaml("config/module_scoring_20261022.yml")
  annotation <- read.csv(rule$annotation_path, stringsAsFactors = FALSE)
  folds <- validate_gse65391_subject_folds_v2(metadata, cohort$subject_folds)
  roles <- read.csv("data/derived/GSE65391/execution_20261022_v1/exploratory_subject_roles.csv",
    stringsAsFactors = FALSE)
  exploratory <- validate_exploratory_subject_roles_v2(metadata, roles)
  expression_path <- "data/interim/GSE65391/expression_processed_probe.rds"
  expression <- NULL
  if (load_expression) {
    expression <- readRDS(expression_path)
    if (anyNA(match(rownames(metadata), colnames(expression))) ||
        anyDuplicated(rownames(expression)) || anyDuplicated(colnames(expression)))
      stop("The recovered expression identifiers do not match the cohort.")
    expression <- expression[, rownames(metadata), drop = FALSE]
    if (any(!is.finite(expression))) stop("The observed expression contains nonfinite values.")
  }
  list(cohort = cohort, metadata = metadata, modules = modules, members = members,
    rule = rule, annotation = annotation, expression = expression, roles = roles,
    fold_audit = folds, exploratory_audit = exploratory,
    source = list(cohort = cohort_path, expression = expression_path,
      recorded_expression = cohort$expression_file, annotation = rule$annotation_path))
}

empirical_reference_base <- function(input, training_rows, filter_size = 10000L) {
  metadata <- input$metadata
  if (is.null(input$expression)) stop("The complete observed expression input is required.")
  design <- if (filter_size > 0L) {
    prepare_gse65391_prediction_design_v2(metadata, input$expression,
      input$annotation, input$modules, input$members, input$rule, training_rows,
      include_sledai = FALSE, filter_size = filter_size)
  } else prepare_gse65391_lowdim_design_v2(metadata, input$expression, input$annotation,
    input$modules, training_rows, input$rule, FALSE, input$members)
  list(clinical = list(A = design$clinical,
    B = prepare_clinical_design_split_v2(metadata, training_rows, include_sledai = TRUE)),
    module_scores = design$X[, grepl("^module__", colnames(design$X)), drop = FALSE],
    genes = design$X[, grepl("^gene__", colnames(design$X)), drop = FALSE],
    training_rows = as.integer(training_rows),
    training_subjects = unique(as.character(metadata$subject_id[training_rows])),
    sample_ids = rownames(metadata), highdim_preprocessing = design$highdim_preprocessing,
    gene_preprocessing = design$gene_preprocessing, source = input$source,
    processing = "reference")
}

empirical_base_path <- function(processing, key) {
  file.path("data/derived/GSE65391/empirical_recovery_20261030/preprocessing",
    processing, key, "prepared_base.rds")
}

empirical_validate_base <- function(base, metadata, training_rows, processing) {
  if (!identical(base$sample_ids, rownames(metadata)) ||
      !identical(base$training_rows, as.integer(training_rows)) ||
      !identical(base$processing, processing) ||
      !setequal(base$training_subjects, metadata$subject_id[training_rows]) ||
      ncol(base$genes) < (if (processing == "reference") 10000L else 5000L) ||
      nrow(base$genes) != nrow(metadata) ||
      ncol(base$module_scores) != 5L || any(!is.finite(base$module_scores)) ||
      any(!is.finite(base$genes))) stop("The saved preprocessing definition differs.")
  invisible(base)
}

empirical_prediction_roles <- function(input, repetition, fold) {
  assignments <- input$cohort$subject_folds
  heldout <- assignments$subject_id[assignments$repetition == repetition & assignments$fold == fold]
  outer <- which(!input$metadata$subject_id %in% heldout)
  inner <- read.csv("data/derived/GSE65391/execution_20261022_v1/inner_subject_folds.csv",
    stringsAsFactors = FALSE)
  inner <- inner[inner$outer_repetition == repetition & inner$outer_fold == fold, ]
  if (!length(heldout) || length(intersect(heldout, inner$subject_id)) ||
      !setequal(inner$subject_id, input$metadata$subject_id[outer]) ||
      anyDuplicated(inner$subject_id) || !setequal(inner$inner_fold, 1:5))
    stop("The inner and outer subject assignments differ from the frozen design.")
  training <- lapply(0:5, function(g) {
    if (g == 0L) outer else which(!input$metadata$subject_id %in%
      c(heldout, inner$subject_id[inner$inner_fold == g]))
  })
  list(heldout = heldout, inner = inner, training = training)
}

empirical_prepare_prediction <- function(input, processing, repetition, fold) {
  roles <- empirical_prediction_roles(input, repetition, fold)
  for (g in 0:5) {
    suffix <- if (g == 0L) "outer" else sprintf("inner%02d", g)
    key <- sprintf("prediction_r%02d_f%02d_%s", repetition, fold, suffix)
    path <- empirical_base_path(processing, key)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    lock <- filelock::lock(paste0(path, ".lock"), timeout = 60000)
    if (is.null(lock)) stop("Unable to acquire the preprocessing lock.")
    if (!file.exists(path)) {
      if (is.null(input$expression)) input <- empirical_inputs()
      base <- empirical_prepare_base(input, roles$training[[g + 1L]], processing)
      empirical_validate_base(base, input$metadata, roles$training[[g + 1L]], processing)
      submission_atomic_rds(base, path, compress = TRUE)
    } else empirical_validate_base(readRDS(path), input$metadata, roles$training[[g + 1L]], processing)
    filelock::unlock(lock)
  }
  invisible(roles)
}

empirical_prepare_base <- function(input, training_rows, processing,
                                   filter_size = if (processing == "reference") 10000L else 5000L) {
  if (processing == "reference") return(empirical_reference_base(input, training_rows, filter_size))
  if (processing == "processed_median_across_mapped_probes")
    return(empirical_median_base(input, training_rows, filter_size))
  if (processing == "official_R1_R2_normexp_training_quantile") {
    Sys.setenv(R_THREADS = "1")
    raw_path <- "data/interim/GSE65391/raw_calibrated_20261030.rds"
    raw <- readRDS(raw_path)
    values <- raw$expression[, rownames(input$metadata), drop = FALSE]
    target <- preprocessCore::normalize.quantiles.determine.target(values[, training_rows, drop = FALSE])
    normalised <- preprocessCore::normalize.quantiles.use.target(values, target)
    dimnames(normalised) <- dimnames(values)
    if (any(!is.finite(normalised)) || any(normalised <= 0))
      stop("The raw-array quantile normalisation produced invalid values.")
    input$expression <- log2(normalised)
    base <- empirical_reference_base(input, training_rows, filter_size)
    base$processing <- processing
    base$source$raw_expression <- raw_path
    base$processing_transformation <- list(target = target,
      target_subjects = unique(as.character(input$metadata$subject_id[training_rows])),
      package = as.character(packageVersion("preprocessCore")))
    return(base)
  }
  stop("Unknown processing route: ", processing)
}

empirical_gene_medians <- function(input) {
  path <- "data/interim/GSE65391/gene_probe_medians_20261030.rds"
  if (!file.exists(path)) {
    filtered <- remove_probe_rows_v2(input$expression, input$annotation)
    groups <- split(seq_len(nrow(filtered$expression)), filtered$annotation$gene_symbol)
    values <- t(vapply(groups, function(rows) {
      if (length(rows) == 1L) as.numeric(filtered$expression[rows, ]) else
        matrixStats::colMedians(filtered$expression[rows, , drop = FALSE])
    }, numeric(ncol(filtered$expression))))
    colnames(values) <- colnames(filtered$expression)
    submission_atomic_rds(list(expression = values, probe_groups = lapply(groups,
      function(rows) rownames(filtered$expression)[rows]), source = input$source), path, compress = TRUE)
  }
  result <- readRDS(path)
  if (anyNA(match(rownames(input$metadata), colnames(result$expression))))
    stop("The median-probe measurements omit required visits.")
  result$expression <- result$expression[, rownames(input$metadata), drop = FALSE]
  result
}

empirical_median_base <- function(input, training_rows, filter_size) {
  median <- empirical_gene_medians(input)
  order <- unlist(input$rule$module_order, use.names = FALSE)
  scores <- matrix(NA_real_, nrow(input$metadata), length(order),
    dimnames = list(rownames(input$metadata), paste0("module__", order)))
  audit <- list()
  for (j in seq_along(order)) {
    members <- input$members[input$members$module_id == order[j], ]
    kind <- unique(members$member_type)
    if (length(kind) != 1L) stop("Module members have inconsistent measurement types.")
    values <- if (kind == "gene_symbol") {
      median$expression[match(members$member_id, rownames(median$expression)), , drop = FALSE]
    } else input$expression[match(members$member_id, rownames(input$expression)), , drop = FALSE]
    if (any(!is.finite(values)) || max(abs(members$weight - 1 / nrow(members))) > 1e-12)
      stop("The frozen module membership is incomplete under median-probe processing.")
    scores[, j] <- as.numeric(crossprod(members$weight, values))
    audit[[j]] <- data.frame(module_id = order[j], member_type = kind,
      source_members = nrow(members), observed_members = nrow(values))
  }
  center <- colMeans(scores[training_rows, , drop = FALSE])
  scale <- apply(scores[training_rows, , drop = FALSE], 2L, sd)
  if (any(!is.finite(scale)) || any(scale <= 0)) stop("The training module variation is insufficient.")
  module_scores <- sweep(sweep(scores, 2L, center, "-"), 2L, scale, "/")
  gene <- NULL
  genes <- matrix(numeric(), nrow(input$metadata), 0L)
  if (filter_size > 0L) {
    selected <- filter_top_training_mad_v2(median$expression, rownames(input$metadata)[training_rows], filter_size)
    gene <- scale_from_training_v2(t(selected$expression), training_rows)
    genes <- gene$X
    colnames(genes) <- paste0("gene__", colnames(genes))
    if (ncol(genes) != filter_size) stop("The median-probe gene filter has insufficient finite variation.")
  }
  list(clinical = list(A = prepare_clinical_design_split_v2(input$metadata, training_rows, FALSE),
    B = prepare_clinical_design_split_v2(input$metadata, training_rows, TRUE)),
    module_scores = module_scores, genes = genes, training_rows = as.integer(training_rows),
    training_subjects = unique(as.character(input$metadata$subject_id[training_rows])),
    sample_ids = rownames(input$metadata), processing = "processed_median_across_mapped_probes",
    source = input$source, gene_preprocessing = list(probe_groups = median$probe_groups),
    highdim_preprocessing = if (is.null(gene)) NULL else
      list(filter_size = filter_size, gene_ids = colnames(genes), center = gene$center, scale = gene$scale),
    module_preprocessing = list(raw_scores = scores, center = center, scale = scale,
      audit = do.call(rbind, audit), score_members = input$members))
}

empirical_processing_module_base <- function(input, processing, omitted = NA_character_, index = 0L) {
  key <- if (index == 0L) "module_full_sample" else sprintf("module_delete_subject_%04d", index)
  path <- empirical_base_path(processing, key)
  keep <- if (is.na(omitted)) rep(TRUE, nrow(input$metadata)) else input$metadata$subject_id != omitted
  expected <- rownames(input$metadata)[keep]
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  lock <- filelock::lock(paste0(path, ".lock"), timeout = 60000)
  if (is.null(lock)) stop("Unable to acquire the module preprocessing lock.")
  on.exit(filelock::unlock(lock), add = TRUE)
  if (!file.exists(path)) {
    input$metadata <- input$metadata[keep, , drop = FALSE]
    input$expression <- input$expression[, keep, drop = FALSE]
    base <- empirical_prepare_base(input, seq_len(nrow(input$metadata)), processing, filter_size = 0L)
    base$omitted_subject <- omitted
    submission_atomic_rds(base, path, compress = TRUE)
  }
  base <- readRDS(path)
  if (!identical(base$omitted_subject, omitted) || !identical(base$processing, processing) ||
      !identical(base$sample_ids, expected) ||
      length(base$training_rows) != length(base$sample_ids) || anyDuplicated(base$sample_ids))
    stop("The module deletion preprocessing differs from the declared role.")
  base
}

empirical_exploratory_base <- function(input, processing, role) {
  audit <- input$exploratory_audit
  selected <- audit[audit$role_id == role, ]
  if (nrow(selected) != 1L) stop("The exploratory role is not registered.")
  assignment <- input$roles[input$roles$split == selected$split &
    input$roles$swapped == selected$swapped, ]
  training <- which(input$metadata$subject_id %in% assignment$subject_id[assignment$role == "selection"])
  inference <- setdiff(seq_len(nrow(input$metadata)), training)
  path <- empirical_base_path(processing, paste0("exploratory_", role))
  if (!file.exists(path)) {
    if (processing == "reference") {
      base <- prepare_exploratory_gene_split_v2(input$metadata, input$expression,
        input$annotation, assignment, filter_size = 5000L)
    } else {
      prepared <- empirical_prepare_base(input, training, processing)
      base <- list(genes = prepared$genes[, 1:5000, drop = FALSE], clinical = prepared$clinical,
        selection_rows = training, inference_rows = inference,
        selection_subjects = unique(as.character(input$metadata$subject_id[training])),
        inference_subjects = unique(as.character(input$metadata$subject_id[inference])),
        sample_ids = rownames(input$metadata), preprocessing = prepared$highdim_preprocessing,
        assignment = assignment)
    }
    base$processing <- processing
    base$source_path <- path
    base$source_inputs <- input$source
    submission_atomic_rds(base, path, compress = FALSE)
  }
  base <- readRDS(path)
  if (!identical(base$sample_ids, rownames(input$metadata)) ||
      !identical(base$selection_rows, training) || !identical(base$inference_rows, inference) ||
      !identical(base$assignment, assignment) || !identical(base$processing, processing) ||
      length(intersect(base$selection_subjects, base$inference_subjects)))
    stop("The saved exploratory preprocessing uses different subjects.")
  base
}
