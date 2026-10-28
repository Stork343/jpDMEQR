assemble_gpl10558_probe_annotation_v2 <- function(probe_ids, platform, vendor, source, build_date) {
  if (!all(c("ID", "Symbol", "Entrez_Gene_ID", "Definition", "Synonyms", "ILMN_Gene", "Species") %in%
           names(vendor)) || nrow(vendor) != nrow(platform) || !setequal(vendor$ID, platform$ID)) {
    stop("Original submitter and curated GPL10558 probe inventories do not reconcile.")
  }
  annotation <- build_gpl10558_probe_annotation_v2(probe_ids, platform, source, build_date)
  position <- match(annotation$probe_id, vendor$ID)
  if (anyNA(position)) stop("Expression contains a probe absent from official GPL10558 sources.")
  annotation$submitter_symbol <- trimws(as.character(vendor$Symbol[position]))
  annotation$submitter_gene_id <- as.character(vendor$Entrez_Gene_ID[position])
  annotation$submitter_gene_name <- as.character(vendor$Definition[position])
  annotation$submitter_synonyms <- as.character(vendor$Synonyms[position])
  annotation$submitter_ilmn_gene <- as.character(vendor$ILMN_Gene[position])
  annotation$gene_symbol_curated <- annotation$gene_symbol
  annotation$is_control <- as.character(vendor$Species[position]) != "Homo sapiens"
  annotation$control_status <- ifelse(annotation$is_control, "submitter_nonhuman_or_no_species", "human_transcript")
  use_submitter <- is.na(annotation$gene_symbol) & !is.na(annotation$submitter_symbol) &
    nzchar(annotation$submitter_symbol) & !annotation$is_control
  annotation$gene_symbol[use_submitter] <- annotation$submitter_symbol[use_submitter]
  annotation$gene_name[use_submitter] <- annotation$submitter_gene_name[use_submitter]
  annotation$entrez_id[use_submitter] <- annotation$submitter_gene_id[use_submitter]
  annotation$annotation_source[use_submitter] <- "GPL10558_original_submitter_table"
  annotation$is_multimapping <- (!is.na(annotation$gene_symbol) &
                                  grepl("///", annotation$gene_symbol, fixed = TRUE)) |
    (!is.na(annotation$entrez_id) & grepl("///", annotation$entrez_id, fixed = TRUE))
  annotation$mapping_status <- ifelse(annotation$is_control, "array_control",
                                     ifelse(is.na(annotation$gene_symbol) | !nzchar(annotation$gene_symbol),
                                             "gene_unannotated",
                                             ifelse(annotation$is_multimapping, "multiple_genes", "unique_gene")))
  annotation
}

read_frozen_module_score_members_v2 <- function(member_path,
                                                freeze_path = paste0(member_path, ".freeze.csv")) {
  if (!file.exists(member_path) || !file.exists(freeze_path)) {
    stop("MODULE_GATE_FAILED: the versioned score members and freeze record are required.")
  }
  members <- utils::read.csv(member_path, stringsAsFactors = FALSE, check.names = FALSE)
  freeze <- utils::read.csv(freeze_path, stringsAsFactors = FALSE)
  required <- c("module_id", "member_type", "member_id", "source_member_id", "weight",
                "source", "source_version", "outcome_association_inspected")
  if (!nrow(members) || !all(required %in% names(members)) || nrow(freeze) != 1L ||
      !all(c("md5", "source", "frozen_at_utc") %in% names(freeze)) ||
      !identical(as.character(unname(tools::md5sum(member_path))), as.character(freeze$md5))) {
    stop("MODULE_GATE_FAILED: the score-member file does not match its complete freeze record.")
  }
  text_fields <- setdiff(required, c("weight", "outcome_association_inspected"))
  if (any(vapply(members[text_fields], function(x) anyNA(x) || any(!nzchar(x)), logical(1))) ||
      anyDuplicated(members[c("module_id", "member_id")]) ||
      anyNA(members$outcome_association_inspected) || any(members$outcome_association_inspected) ||
      any(!members$member_type %in% c("gene_symbol", "ILMN_probe")) ||
      any(!is.finite(members$weight)) || any(members$weight <= 0)) {
    stop("MODULE_GATE_FAILED: invalid source-faithful score membership.")
  }
  for (id in unique(members$module_id)) {
    rows <- members[members$module_id == id, ]
    if (length(unique(rows$member_type)) != 1L ||
        max(abs(rows$weight - 1 / nrow(rows))) > 1e-12) {
      stop("MODULE_GATE_FAILED: each module requires one representation and complete equal weights.")
    }
  }
  members
}

score_frozen_module_representations_v2 <- function(expression_probe, annotation,
                                                   module_signatures, score_members,
                                                   training_rows, module_rule) {
  expression_probe <- assert_numeric_matrix_v2(expression_probe, "expression_probe")
  annotation <- as.data.frame(annotation, stringsAsFactors = FALSE)
  module_signatures <- as.data.frame(module_signatures, stringsAsFactors = FALSE)
  score_members <- as.data.frame(score_members, stringsAsFactors = FALSE)
  required_rule <- c("module_order", "representations", "minimum_fraction", "score",
                     "score_transformation", "probe_collapse_rule")
  if (!all(required_rule %in% names(module_rule)) || module_rule$minimum_fraction != 1 ||
      module_rule$score != "equal_weight_mean_processed_normalised_expression" ||
      module_rule$score_transformation != "module_mean_center_and_sd_from_training_rows" ||
      module_rule$probe_collapse_rule != "highest_training_mad") {
    stop("MODULE_RULE_REQUIRED: the complete frozen representation and scoring rule is required.")
  }
  sample_ids <- colnames(expression_probe)
  if (is.null(sample_ids) || anyDuplicated(sample_ids) || anyNA(sample_ids) ||
      is.null(rownames(expression_probe)) || anyDuplicated(rownames(expression_probe))) {
    stop("MODULE_GATE_FAILED: expression sample and probe identifiers must be unique.")
  }
  train <- if (is.logical(training_rows)) {
    if (length(training_rows) != ncol(expression_probe) || anyNA(training_rows)) {
      stop("Invalid training_rows for module scoring.")
    }
    which(training_rows)
  } else as.integer(training_rows)
  if (length(train) < 2L || anyNA(train) || any(train < 1L | train > ncol(expression_probe)) ||
      anyDuplicated(train)) stop("Module scoring requires at least two distinct training rows.")
  order <- as.character(unlist(module_rule$module_order))
  if (!length(order) || anyDuplicated(order) || !setequal(order, score_members$module_id) ||
      !all(c("probe_id", "gene_symbol", "is_control", "is_multimapping") %in% names(annotation)) ||
      anyDuplicated(annotation$probe_id)) {
    stop("MODULE_GATE_FAILED: score members, rule order and annotation are inconsistent.")
  }
  genes <- unique(score_members$member_id[score_members$member_type == "gene_symbol"])
  relevant <- !is.na(annotation$gene_symbol) & annotation$gene_symbol %in% genes &
    !is.na(annotation$is_control) & !annotation$is_control &
    !is.na(annotation$is_multimapping) & !annotation$is_multimapping
  positions <- match(annotation$probe_id[relevant], rownames(expression_probe))
  present <- !is.na(positions)
  gene_annotation <- annotation[relevant, , drop = FALSE][present, , drop = FALSE]
  gene_expression <- expression_probe[positions[present], , drop = FALSE]
  if (!setequal(genes, unique(gene_annotation$gene_symbol)) || any(!is.finite(gene_expression))) {
    stop("MODULE_GATE_FAILED: every gene member requires observed finite expression.")
  }
  collapsed <- collapse_probes_highest_training_mad_v2(
    gene_expression, gene_annotation$gene_symbol, sample_ids[train]
  )
  gene_probe_map <- data.frame(
    gene_symbol = rownames(collapsed$expression),
    probe_id = rownames(gene_expression)[collapsed$selected_probe_rows],
    training_mad = collapsed$training_mad, stringsAsFactors = FALSE
  )
  raw_scores <- matrix(NA_real_, ncol(expression_probe), length(order),
                        dimnames = list(sample_ids, order))
  audit <- vector("list", length(order))
  for (i in seq_along(order)) {
    id <- order[i]
    rows <- score_members[score_members$module_id == id, , drop = FALSE]
    unit <- unique(rows$member_type)
    if (length(unit) != 1L || !identical(unit, as.character(module_rule$representations[[id]])) ||
        anyDuplicated(rows$member_id) || any(!is.finite(rows$weight)) ||
        max(abs(rows$weight - 1 / nrow(rows))) > 1e-12) {
      stop("MODULE_GATE_FAILED: module representation or member weights differ from the frozen rule.")
    }
    if (unit == "gene_symbol") {
      source <- module_signatures[module_signatures$module_id == id, , drop = FALSE]
      if (!setequal(rows$member_id, source$gene_symbol) || nrow(rows) != nrow(source)) {
        stop("MODULE_GATE_FAILED: the gene representation must retain all public members.")
      }
      values <- collapsed$expression[match(rows$member_id, rownames(collapsed$expression)), , drop = FALSE]
    } else if (unit == "ILMN_probe") {
      index <- match(rows$member_id, rownames(expression_probe))
      if (anyNA(index)) stop("MODULE_GATE_FAILED: an original module probe is absent from expression.")
      values <- expression_probe[index, , drop = FALSE]
    } else stop("MODULE_GATE_FAILED: unknown module measurement unit.")
    if (any(!is.finite(values))) stop("MODULE_GATE_FAILED: module expression contains missing measurements.")
    raw_scores[, i] <- as.numeric(crossprod(rows$weight, values))
    audit[[i]] <- data.frame(module_id = id, member_type = unit,
                            source_members = nrow(rows), observed_members = nrow(values),
                            fraction_observed = 1, weight_per_member = rows$weight[1],
                            eligible = TRUE, stringsAsFactors = FALSE)
  }
  center <- colMeans(raw_scores[train, , drop = FALSE])
  scale <- apply(raw_scores[train, , drop = FALSE], 2, stats::sd)
  if (any(!is.finite(scale)) || any(scale <= 0)) {
    stop("MODULE_GATE_FAILED: every module must have positive training expression variation.")
  }
  scores <- sweep(sweep(raw_scores, 2, center, "-"), 2, scale, "/")
  list(scores = scores, raw_scores = raw_scores, audit = do.call(rbind, audit),
       center = center, scale = scale, gene_probe_map = gene_probe_map,
       score_members = score_members, training_sample_ids = sample_ids[train])
}
