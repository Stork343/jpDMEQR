audit_gse65391_identifiers_v2 <- function(metadata, expression_sample_ids = NULL) {
  required <- c("geo_accession", "title", "subject_id", "visit_number", "visit_id", "cumulative_time",
                "group_type", "technical_replicate_flag", "batch_replicate_flag")
  missing <- setdiff(required, names(metadata))
  if (length(missing)) stop("Metadata fields missing: ", paste(missing, collapse = ", "))
  md <- as.data.frame(metadata, stringsAsFactors = FALSE)
  records <- list()
  add_conflict <- function(rows, issue) {
    if (length(rows)) records[[length(records) + 1L]] <<- data.frame(
      sample_id = as.character(md$geo_accession[rows]),
      issue = issue, stringsAsFactors = FALSE
    )
  }
  duplicate_all <- function(x) duplicated(x) | duplicated(x, fromLast = TRUE)
  biological <- md$group_type == "SLE" & !md$technical_replicate_flag
  biological[is.na(biological)] <- FALSE
  sample_id <- standardize_geo_missing_v2(md$geo_accession)
  subject_id <- standardize_geo_missing_v2(md$subject_id)
  visit <- suppressWarnings(as.numeric(md$visit_number))
  add_conflict(which(is.na(sample_id)), "missing_sample_id")
  add_conflict(which(duplicate_all(sample_id)), "duplicate_sample_id")
  add_conflict(which(is.na(md$group_type) | md$group_type == "UNRESOLVED"), "unresolved_cohort")
  add_conflict(which(is.na(md$technical_replicate_flag) | is.na(md$batch_replicate_flag)),
               "unresolved_replicate_status")
  add_conflict(which(biological & is.na(subject_id)), "missing_subject_id")
  valid_visit <- is.finite(visit) & visit >= 1 & visit == floor(visit)
  valid_time <- is.finite(md$cumulative_time)
  subject_time <- paste(subject_id, sprintf("%.17g", md$cumulative_time), sep = "::")
  time_unique <- !duplicate_all(subject_time)
  sample_unique <- !is.na(sample_id) & !duplicate_all(sample_id)
  biological_id <- paste(subject_id, sprintf("V%03d", visit), sep = "_")
  add_conflict(which(biological & valid_visit & md$visit_id != biological_id), "visit_id_mismatch")
  bio_rows <- which(biological)

  title_subject <- sub(".*-(SLE-[0-9]+)-V.*", "\\1", md$title)
  title_visit <- suppressWarnings(as.numeric(sub(".*-V([0-9]+)-.*", "\\1", md$title)))
  subject_parsed <- grepl("-(SLE-[0-9]+)-V", md$title)
  add_conflict(which(biological & !subject_parsed), "title_subject_unresolved")
  time_key_allowed <- biological & !valid_visit & valid_time & time_unique &
    sample_unique & !is.na(subject_id) & subject_parsed & title_subject == subject_id
  add_conflict(which(biological & !valid_visit & !time_key_allowed), "unresolved_visit_identity")
  add_conflict(which(biological & valid_visit & !is.finite(title_visit)), "title_visit_unresolved")
  add_conflict(which(biological & subject_parsed & !is.na(subject_id) &
                       title_subject != subject_id), "title_subject_mismatch")
  add_conflict(which(biological & is.finite(title_visit) & is.finite(visit) &
                       title_visit != visit), "title_visit_mismatch")
  key_source <- ifelse(valid_visit, "subject_visit", ifelse(time_key_allowed, "subject_time", "unresolved"))
  canonical_key <- ifelse(key_source == "subject_visit", paste0("subject_visit::", biological_id),
                         ifelse(key_source == "subject_time", paste0("subject_time::", subject_time), NA_character_))
  add_conflict(bio_rows[duplicate_all(canonical_key[bio_rows])], "duplicate_biological_visit")
  canonical_keys <- data.frame(
    sample_id = sample_id, biological = biological, canonical_visit_key = canonical_key,
    key_source = key_source, recorded_visit = visit, recorded_cumulative_time = md$cumulative_time,
    stringsAsFactors = FALSE
  )
  expression_checked <- !is.null(expression_sample_ids)
  expression_problems <- character()
  if (expression_checked) {
    expression_sample_ids <- as.character(expression_sample_ids)
    if (anyNA(expression_sample_ids) || anyDuplicated(expression_sample_ids)) {
      expression_problems <- c(expression_problems, "invalid_expression_sample_ids")
    }
    if (!setequal(expression_sample_ids, sample_id) ||
        length(expression_sample_ids) != nrow(md)) {
      expression_problems <- c(expression_problems, "expression_metadata_mismatch")
    }
  }
  conflicts <- if (length(records)) do.call(rbind, records) else data.frame(
    sample_id = character(), issue = character(), stringsAsFactors = FALSE
  )
  issue_counts <- if (nrow(conflicts)) as.data.frame(table(conflicts$issue)) else
    data.frame(Var1 = character(), Freq = integer())
  names(issue_counts) <- c("issue", "sample_count")
  list(
    conflicts = conflicts, issue_counts = issue_counts, canonical_keys = canonical_keys,
    expression_problems = expression_problems,
    summary = data.frame(
      biological_visits = sum(biological),
      flagged_sle_batch_records = sum(md$group_type == "SLE" & md$batch_replicate_flag),
      subject_time_keys = sum(biological & key_source == "subject_time"),
      critical_conflicts = nrow(conflicts) + length(expression_problems),
      expression_identity_checked = expression_checked,
      pass = nrow(conflicts) == 0L && !length(expression_problems) && expression_checked
    )
  )
}

audit_gse65391_replicate_relations_v2 <- function(metadata, source_metadata = NULL) {
  md <- as.data.frame(metadata, stringsAsFactors = FALSE)
  flagged <- which(md$batch_replicate_flag | md$technical_replicate_flag)
  rows <- lapply(flagged, function(ii) {
    other <- setdiff(seq_len(nrow(md)), ii)
    same_subject <- !is.na(md$subject_id[other]) & !is.na(md$subject_id[ii]) &
      md$subject_id[other] == md$subject_id[ii]
    same_visit <- same_subject & is.finite(md$visit_number[other]) &
      is.finite(md$visit_number[ii]) & md$visit_number[other] == md$visit_number[ii]
    same_time <- same_subject & is.finite(md$cumulative_time[other]) &
      is.finite(md$cumulative_time[ii]) & md$cumulative_time[other] == md$cumulative_time[ii]
    source_flag <- NA_character_
    physical_matches <- NA_integer_
    if (!is.null(source_metadata)) {
      source_row <- match(md$geo_accession[ii], rownames(source_metadata))
      if (is.na(source_row)) stop("Original GEO phenotype row is absent.")
      if ("batch_replicate:ch1" %in% names(source_metadata)) {
        source_flag <- as.character(source_metadata[source_row, "batch_replicate:ch1"])
      }
      if ("description" %in% names(source_metadata)) {
        physical_matches <- sum(source_metadata$description ==
                                  source_metadata$description[source_row], na.rm = TRUE) - 1L
      }
    }
    data.frame(
      sample_id = md$geo_accession[ii], group_type = md$group_type[ii],
      batch_replicate_flag = md$batch_replicate_flag[ii],
      technical_replicate_flag = md$technical_replicate_flag[ii],
      original_batch_replicate = source_flag,
      matching_subject_visit_records = sum(same_visit, na.rm = TRUE),
      matching_subject_time_records = sum(same_time, na.rm = TRUE),
      matching_physical_array_records = physical_matches,
      relation_status = if (sum(same_visit, na.rm = TRUE) > 0L)
        "matched_biological_visit" else if (md$group_type[ii] == "SLE")
          "sole_observed_biological_record" else "unmatched_flagged_record",
      stringsAsFactors = FALSE
    )
  })
  relations <- if (length(rows)) do.call(rbind, rows) else data.frame()
  counts <- if (nrow(relations)) stats::aggregate(
    list(records = rep(1L, nrow(relations))),
    relations[c("group_type", "batch_replicate_flag", "technical_replicate_flag", "relation_status")],
    sum
  ) else data.frame()
  list(records = relations, summary = counts,
       source_processing = if (is.null(source_metadata)) character() else
         unique(as.character(source_metadata$data_processing)))
}

gse65391_exclusion_reasons_v2 <- function(metadata, eligible_visit) {
  reason <- rep("included", nrow(metadata))
  assign_reason <- function(condition, label) {
    ii <- which(condition & reason == "included")
    reason[ii] <<- label
  }
  assign_reason(is.na(metadata$group_type) | metadata$group_type == "UNRESOLVED", "unresolved_cohort")
  assign_reason(metadata$group_type == "HEALTHY" & metadata$technical_replicate_flag,
                "healthy_technical_replicate")
  assign_reason(metadata$group_type == "HEALTHY", "healthy_biological_sample")
  assign_reason(metadata$technical_replicate_flag, "technical_replicate")
  assign_reason(!is.finite(metadata$c3), "missing_c3")
  assign_reason(is.finite(metadata$c3) & metadata$c3 <= 0, "nonpositive_c3")
  assign_reason(is.na(standardize_geo_missing_v2(metadata$subject_id)), "missing_subject_id")
  assign_reason(!eligible_visit, "fewer_than_two_eligible_visits")
  data.frame(sample_id = metadata$geo_accession, reason = reason, stringsAsFactors = FALSE)
}

audit_gse65391_batch_outcome_v2 <- function(metadata, eligible_visit, quantiles) {
  md <- metadata[eligible_visit, , drop = FALSE]
  y <- log(as.numeric(md$c3))
  if (!length(y)) stop("No eligible C3 values for the batch audit.")
  breaks <- unique(as.numeric(stats::quantile(y, quantiles, names = FALSE)))
  region <- cut(y, breaks = c(-Inf, breaks, Inf), include.lowest = TRUE)
  batch <- standardize_geo_missing_v2(md$batch)
  missing_batch <- is.na(batch)
  batch[missing_batch] <- "Missing"
  counts <- as.data.frame(table(batch = batch, outcome_region = region))
  summary <- do.call(rbind, lapply(unique(batch), function(b) {
    ix <- batch == b
    data.frame(batch = b, visits = sum(ix), subjects = length(unique(md$subject_id[ix])),
               outcome_regions = length(unique(region[ix])),
               outcome_min = min(y[ix]), outcome_max = max(y[ix]),
               pass = b != "Missing" && length(unique(region[ix])) > 1L,
               stringsAsFactors = FALSE)
  }))
  list(counts = counts, summary = summary, pass = all(summary$pass))
}

bootstrap_gse65391_curvature_v2 <- function(metadata, diagnostics, rule) {
  required <- c("version", "method", "source", "primary_quantiles", "sensitivity_quantiles",
                "h_multipliers", "family_confidence", "bootstrap_replicates", "bootstrap_seed",
                "quantile_type", "dimension_rule", "minimum_subject_rule", "sampling_unit",
                "recompute_quantiles", "cluster_weighting")
  if (!all(required %in% names(rule)) || !identical(rule$method, "subject_bootstrap_percentile") ||
      !identical(rule$sampling_unit, "subject") || !isTRUE(rule$recompute_quantiles) ||
      !identical(rule$cluster_weighting, "equal_subject") ||
      !identical(rule$dimension_rule, "floor_n_over_log_n") ||
      !identical(rule$minimum_subject_rule, "d_max_plus_one")) {
    stop("Invalid frozen subject-bootstrap curvature rule.")
  }
  B <- as.integer(rule$bootstrap_replicates)
  confidence <- as.numeric(rule$family_confidence)
  if (B != 10000L || rule$bootstrap_seed != 20261022L || confidence != 0.95 ||
      rule$quantile_type != 7L || !nzchar(rule$source) ||
      !setequal(unique(diagnostics$tau[diagnostics$primary]), unlist(rule$primary_quantiles)) ||
      !setequal(unique(diagnostics$tau[!diagnostics$primary]), unlist(rule$sensitivity_quantiles)) ||
      !setequal(unique(diagnostics$h_multiplier), unlist(rule$h_multipliers))) {
    stop("Curvature diagnostics differ from the versioned bootstrap rule.")
  }
  y <- log(as.numeric(metadata$c3))
  subject <- match(as.character(metadata$subject_id), unique(as.character(metadata$subject_id)))
  n <- max(subject)
  visit_count <- tabulate(subject, nbins = n)
  row_weight <- 1 / visit_count[subject]
  taus <- unique(diagnostics$tau)
  tau_index <- match(diagnostics$tau, taus)
  means <- matrix(NA_real_, B, nrow(diagnostics))
  set.seed(as.integer(rule$bootstrap_seed))
  for (b in seq_len(B)) {
    selected <- tabulate(sample.int(n, n, replace = TRUE), nbins = n)
    boot_y <- rep(y, selected[subject])
    locations <- as.numeric(stats::quantile(boot_y, taus, type = rule$quantile_type, names = FALSE))
    boot_weight <- selected[subject] * row_weight / n
    for (j in seq_len(nrow(diagnostics))) {
      means[b, j] <- sum(boot_weight * phi_smooth_v2(
        y - locations[tau_index[j]], diagnostics$tau[j], diagnostics$h[j]
      ))
    }
  }
  diagnostics$family <- ifelse(diagnostics$primary, "primary", "sensitivity")
  diagnostics$family_size <- ifelse(diagnostics$primary, sum(diagnostics$primary),
                                    sum(!diagnostics$primary))
  diagnostics$lower_percentile <- (1 - confidence) / diagnostics$family_size
  diagnostics$bootstrap_lower_bound <- vapply(seq_len(nrow(diagnostics)), function(j) {
    as.numeric(stats::quantile(means[, j], diagnostics$lower_percentile[j],
                              type = rule$quantile_type, names = FALSE))
  }, numeric(1))
  d_max <- floor(n / log(max(n, 3)))
  diagnostics$d_max <- d_max
  diagnostics$minimum_subjects <- d_max + 1L
  diagnostics$bootstrap_replicates <- B
  diagnostics$bootstrap_seed <- as.integer(rule$bootstrap_seed)
  diagnostics$family_confidence <- confidence
  diagnostics$pass <- diagnostics$bootstrap_lower_bound > 0 &
    diagnostics$local_subjects >= diagnostics$minimum_subjects &
    diagnostics$local_subject_effective_n >= diagnostics$minimum_subjects
  diagnostics$status <- ifelse(diagnostics$pass, "pass", "failed")
  list(diagnostics = diagnostics, bootstrap_means = means, rng_kind = RNGkind(),
       primary_pass = all(diagnostics$pass[diagnostics$primary]),
       sensitivity_pass = all(diagnostics$pass[!diagnostics$primary]))
}

audit_gse65391_local_curvature_v2 <- function(metadata, eligible_visit, quantiles,
                                             primary_quantiles, multipliers,
                                             curvature_rule = NULL) {
  md <- metadata[eligible_visit, , drop = FALSE]
  y <- log(as.numeric(md$c3))
  subjects <- as.character(md$subject_id)
  n <- length(unique(subjects))
  if (!length(y) || n < 1L) stop("No eligible C3 values for the curvature audit.")
  grid <- expand.grid(tau = quantiles, h_multiplier = multipliers)
  diagnostics <- do.call(rbind, lapply(seq_len(nrow(grid)), function(ii) {
    tau <- grid$tau[ii]
    h <- grid$h_multiplier[ii] * n^(-3 / 10)
    location <- as.numeric(stats::quantile(y, tau, names = FALSE))
    weights <- phi_smooth_v2(y - location, tau, h)
    subject_weight <- tapply(weights, subjects, mean)
    data.frame(
      tau = tau, primary = tau %in% primary_quantiles,
      h_multiplier = grid$h_multiplier[ii], h = h, log_c3_quantile = location,
      local_visits = sum(weights > 0),
      local_subjects = length(unique(subjects[weights > 0])),
      local_subject_effective_n = if (sum(subject_weight^2) > 0)
        sum(subject_weight)^2 / sum(subject_weight^2) else 0,
      cluster_average_curvature = mean(subject_weight), stringsAsFactors = FALSE
    )
  }))
  # 缺少冻结规则时保留待审阅状态。
  frozen <- !is.null(curvature_rule)
  if (frozen && identical(curvature_rule$method, "subject_bootstrap_percentile")) {
    bootstrap <- bootstrap_gse65391_curvature_v2(md, diagnostics, curvature_rule)
    return(list(diagnostics = bootstrap$diagnostics, rule = curvature_rule,
                pass = bootstrap$primary_pass, sensitivity_pass = bootstrap$sensitivity_pass,
                bootstrap_means = bootstrap$bootstrap_means, rng_kind = bootstrap$rng_kind))
  }
  if (frozen) {
    required <- c("minimum_local_subjects", "minimum_effective_subjects", "source")
    if (!all(required %in% names(curvature_rule)) ||
        !nzchar(curvature_rule$source) ||
        any(!is.finite(c(curvature_rule$minimum_local_subjects,
                         curvature_rule$minimum_effective_subjects)))) {
      stop("Invalid frozen curvature rule.")
    }
    diagnostics$pass <- diagnostics$local_subjects >= curvature_rule$minimum_local_subjects &
      diagnostics$local_subject_effective_n >= curvature_rule$minimum_effective_subjects &
      diagnostics$cluster_average_curvature > 0
    diagnostics$status <- ifelse(diagnostics$pass, "pass", "failed")
  } else {
    diagnostics$pass <- NA
    diagnostics$status <- "review_required"
  }
  list(diagnostics = diagnostics, rule = curvature_rule,
       pass = frozen && all(diagnostics$pass[diagnostics$primary]))
}

complete_gse65391_audit_v2 <- function(audit, metadata, expression_sample_ids,
                                      primary_quantiles, sensitivity_quantiles,
                                      bandwidth_multipliers, curvature_rule) {
  identifiers <- audit_gse65391_identifiers_v2(metadata, expression_sample_ids)
  batch <- audit_gse65391_batch_outcome_v2(metadata, audit$eligible_visit, primary_quantiles)
  curvature <- audit_gse65391_local_curvature_v2(
    metadata, audit$eligible_visit, unique(c(primary_quantiles, sensitivity_quantiles)),
    primary_quantiles, bandwidth_multipliers, curvature_rule
  )
  added <- data.frame(
    gate = c("subject_visit_reconciliation", "batch_outcome_overlap", "local_curvature"),
    observed = c(identifiers$summary$critical_conflicts,
                 min(batch$summary$outcome_regions),
                 min(curvature$diagnostics$local_subject_effective_n[
                   curvature$diagnostics$primary])),
    threshold = c(0, 1, NA_real_), direction = c("==", ">", "frozen_rule"),
    pass = c(if (identifiers$summary$critical_conflicts > 0L) FALSE else
               if (identifiers$summary$expression_identity_checked)
                 identifiers$summary$pass else NA, batch$pass,
             if (is.null(curvature_rule)) NA else curvature$pass),
    stringsAsFactors = FALSE
  )
  audit$gates$status <- ifelse(audit$gates$pass, "pass", "failed")
  added$status <- ifelse(is.na(added$pass), "review_required",
                         ifelse(added$pass, "pass", "failed"))
  audit$gates <- rbind(audit$gates, added)
  audit$summary$scalar_gates_pass <- isTRUE(audit$summary$all_gates_pass)
  audit$summary$all_gates_pass <- all(!is.na(audit$gates$pass) & audit$gates$pass)
  audit$summary$audit_complete <- !anyNA(audit$gates$pass) &&
    identifiers$summary$expression_identity_checked
  audit$summary$review_required <- anyNA(audit$gates$pass) ||
    !identifiers$summary$expression_identity_checked
  audit$summary$flagged_sle_batch_records <- identifiers$summary$flagged_sle_batch_records
  audit$summary$subject_time_keys <- identifiers$summary$subject_time_keys
  audit$identifiers <- identifiers
  audit$batch <- batch
  audit$curvature <- curvature
  audit$audit_version <- "gse65391_full_v1"
  audit
}

validate_gse65391_fit_inputs_v2 <- function(audit, module_path,
                                            module_freeze_path = paste0(module_path, ".freeze.csv")) {
  required <- c("minimum_subjects", "minimum_visits", "minimum_repeated_subjects",
                "maximum_missing_fraction", "positive_values", "minimum_distinct_values",
                "maximum_modal_fraction", "subject_visit_reconciliation", "batch_outcome_overlap",
                "local_curvature")
  if (!identical(audit$audit_version, "gse65391_full_v1") ||
      !all(required %in% audit$gates$gate) ||
      !isTRUE(audit$summary$audit_complete) ||
      !isTRUE(audit$summary$all_gates_pass)) {
    stop("EMPIRICAL_GATE_INCOMPLETE: complete the frozen C3 audit before fitting.")
  }
  read_frozen_module_definitions_v2(module_path, module_freeze_path)
}

read_frozen_module_definitions_v2 <- function(module_path,
                                               module_freeze_path = paste0(module_path, ".freeze.csv")) {
  if (!file.exists(module_path)) stop("MODULE_GATE_FAILED: module_signatures.csv is missing.")
  modules <- utils::read.csv(module_path, stringsAsFactors = FALSE, check.names = FALSE)
  fields <- c("module_id", "source", "source_version", "gene_symbol", "inclusion_reason")
  if (!all(fields %in% names(modules)) || !nrow(modules)) {
    stop("MODULE_GATE_FAILED: module definitions must be nonempty and contain all required fields.")
  }
  if (any(vapply(modules[fields], function(x) anyNA(x) || any(!nzchar(trimws(x))), logical(1))) ||
      anyDuplicated(modules[c("module_id", "gene_symbol")])) {
    stop("MODULE_GATE_FAILED: module definitions contain missing fields or duplicate memberships.")
  }
  if (!file.exists(module_freeze_path)) stop("MODULE_GATE_FAILED: module freeze record is missing.")
  freeze <- utils::read.csv(module_freeze_path, stringsAsFactors = FALSE)
  if (!all(c("md5", "source", "frozen_at_utc") %in% names(freeze)) || nrow(freeze) != 1L ||
      anyNA(freeze[c("md5", "source", "frozen_at_utc")]) ||
      any(!nzchar(unlist(freeze[c("md5", "source", "frozen_at_utc")])))) {
    stop("MODULE_GATE_FAILED: module freeze record is incomplete.")
  }
  if (!identical(as.character(unname(tools::md5sum(module_path))), as.character(freeze$md5))) {
    stop("MODULE_GATE_FAILED: module definitions differ from the frozen record.")
  }
  modules
}

validate_application_implementation_v2 <- function(root, checks_path, inference = FALSE,
                                                   method_version = NULL) {
  if (!file.exists(checks_path)) {
    stop("METHOD_GATE_INCOMPLETE: unit/package check record is missing.")
  }
  checks <- jsonlite::read_json(checks_path, simplifyVector = TRUE)
  required <- c("unit_checks_pass", "package_checks_pass", "commit_sha")
  if (!all(required %in% names(checks)) || !isTRUE(checks$unit_checks_pass) ||
      !isTRUE(checks$package_checks_pass)) {
    stop("METHOD_GATE_INCOMPLETE: unit/package checks have not passed.")
  }
  if (isTRUE(inference) && !isTRUE(checks$generic_target_review_pass)) {
    stop("TARGET_REVIEW_REQUIRED: the application inference route is awaiting theory review.")
  }
  paired <- identical(method_version, "profile_pair_20261022")
  if (paired && (!identical(checks$method_version, method_version) ||
                 !identical(checks$score_definition, "profile_bandwidth_pair_2_half_minus_full") ||
                 is.null(checks$geometry_manifest) || is.null(checks$micro_manifest))) {
    stop("METHOD_GATE_INCOMPLETE: matching paired-method checks and manifest paths are required.")
  }
  geometry_path <- if (paired) file.path(root, checks$geometry_manifest) else
    file.path(root, "results", "preflight", "geometry", "manifest.json")
  micro_path <- if (paired) file.path(root, checks$micro_manifest) else
    file.path(root, "results", "preflight", "micro_preflight_manifest.json")
  geometry <- jsonlite::read_json(geometry_path,
                                 simplifyVector = TRUE)
  micro <- jsonlite::read_json(micro_path,
                              simplifyVector = TRUE)
  if (paired && (!identical(geometry$method_version, method_version) ||
                 !identical(micro$method_version, method_version) ||
                 !identical(geometry$score_definition, checks$score_definition) ||
                 !identical(micro$score_definition, checks$score_definition))) {
    stop("METHOD_GATE_INCOMPLETE: numerical manifests use a different score definition.")
  }
  commit <- system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"), stdout = TRUE)
  if (length(commit) != 1L || !isTRUE(geometry$pass) || !isTRUE(geometry$strict) ||
      !isTRUE(geometry$dantzig_checked) || !isTRUE(geometry$dantzig_pass) ||
      !isTRUE(micro$pass) || !identical(as.character(geometry$commit_sha), commit) ||
      !identical(as.character(micro$commit_sha), commit) ||
      !identical(as.character(checks$commit_sha), commit)) {
    stop("METHOD_GATE_INCOMPLETE: implementation checks do not certify the current commit.")
  }
  list(commit_sha = commit, checks_path = normalizePath(checks_path),
       method_version = checks$method_version,
       score_definition = checks$score_definition)
}

build_gpl10558_probe_annotation_v2 <- function(probe_ids, platform, source, build_date) {
  required <- c("ID", "Gene symbol", "Gene title", "Gene ID")
  if (!all(required %in% names(platform)) || anyDuplicated(platform$ID) ||
      anyDuplicated(probe_ids) || !nzchar(source) || !nzchar(build_date)) {
    stop("Invalid official GPL10558 annotation input.")
  }
  position <- match(probe_ids, platform$ID)
  symbol <- trimws(as.character(platform[["Gene symbol"]][position]))
  gene_id <- trimws(as.character(platform[["Gene ID"]][position]))
  symbol[is.na(symbol) | symbol %in% c("", "NA", "---")] <- NA_character_
  multiple <- (!is.na(symbol) & grepl("///", symbol, fixed = TRUE)) |
    (!is.na(gene_id) & grepl("///", gene_id, fixed = TRUE))
  data.frame(
    probe_id = as.character(probe_ids), gene_symbol = symbol,
    gene_name = as.character(platform[["Gene title"]][position]), entrez_id = gene_id,
    is_control = FALSE, is_multimapping = multiple,
    mapping_status = ifelse(is.na(position), "probe_absent_from_GPL_annotation",
                            ifelse(is.na(symbol), "gene_unannotated",
                                    ifelse(multiple, "multiple_genes", "unique_gene"))),
    control_status = "no_control_indicator_in_GEO_annotation",
    annotation_source = source, annotation_build_date = build_date,
    stringsAsFactors = FALSE
  )
}
