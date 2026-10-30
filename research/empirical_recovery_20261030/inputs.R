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
  design <- prepare_gse65391_prediction_design_v2(metadata, input$expression,
    input$annotation, input$modules, input$members, input$rule, training_rows,
    include_sledai = FALSE, filter_size = filter_size)
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
      ncol(base$genes) < 10000L || nrow(base$genes) != nrow(metadata) ||
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

empirical_prepare_base <- function(input, training_rows, processing) {
  if (processing == "reference") return(empirical_reference_base(input, training_rows))
  stop("The declared processing input has not yet been prepared: ", processing)
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
