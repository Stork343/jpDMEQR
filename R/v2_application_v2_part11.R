subset_expression_to_frozen_modules_v2 <- function(expression_probe, annotation, module_score_members) {
  genes <- module_score_members$member_id[module_score_members$member_type == "gene_symbol"]
  annotation_rows <- !is.na(annotation$gene_symbol) & annotation$gene_symbol %in% genes &
    !is.na(annotation$is_control) & !annotation$is_control &
    !is.na(annotation$is_multimapping) & !annotation$is_multimapping
  probes <- unique(c(annotation$probe_id[annotation_rows],
                      module_score_members$member_id[module_score_members$member_type == "ILMN_probe"]))
  expression_probe[rownames(expression_probe) %in% probes, , drop = FALSE]
}

prepare_application_matched_population_v2 <- function(
    metadata, expression_probe, annotation, module_signatures, module_score_members,
    module_rule, config, model = "A", q = 1L) {
  if (!model %in% c("A", "B") || !q %in% c(1L, 2L) ||
      length(unique(metadata$subject_id)) != 129L || nrow(metadata) != 867L) {
    stop("The matched population requires the complete declared SLE design.")
  }
  design <- prepare_gse65391_lowdim_design_v2(
    metadata, expression_probe, annotation, module_signatures, rep(TRUE, nrow(metadata)),
    module_rule, include_sledai = model == "B", module_score_members = module_score_members
  )
  X <- design$X
  beta <- setNames(rep(config$generating_coefficients$categorical_coefficients, ncol(X)), colnames(X))
  beta["(Intercept)"] <- config$generating_coefficients$intercept
  for (coordinate in intersect(c("centered_time", "age", "sledai"), colnames(X))) {
    value <- config$generating_coefficients[[paste0(coordinate, "_per_population_sd")]]
    beta[coordinate] <- value / stats::sd(X[, coordinate])
  }
  for (coordinate in grep("^module__", colnames(X), value = TRUE)) {
    beta[coordinate] <- config$generating_coefficients[[coordinate]]
  }
  Z <- if (q == 2L) cbind(intercept = 1, centered_time = X[, "centered_time"]) else design$Z
  random_effect_sd <- config$random_effects$intercept_sd
  if (q == 2L) random_effect_sd <- c(random_effect_sd,
                                      config$random_effects$slope_sd_per_population_time_sd /
                                        stats::sd(X[, "centered_time"]))
  list(X = X, Z = Z, beta0 = beta, model = model, q = q,
       metadata = metadata, rows = split(seq_len(nrow(metadata)), as.character(metadata$subject_id)),
       random_effect_sd = random_effect_sd,
       preprocessing = list(clinical = design$clinical, gene_preprocessing = design$gene_preprocessing,
                              module_audit = design$module_audit),
       population_definition = "equal_subject_empirical_distribution_with_fixed_preprocessing")
}

generate_application_matched_population_data_v2 <- function(population, config, error_family,
                                                            n_subjects, seed) {
  if (!error_family %in% c("normal", "t3")) stop("Unknown declared matched error family.")
  set.seed(seed)
  types <- sample.int(length(population$rows), n_subjects, replace = TRUE)
  rows <- unlist(population$rows[types], use.names = FALSE)
  sizes <- lengths(population$rows[types])
  cluster <- rep(seq_len(n_subjects), sizes)
  b <- matrix(stats::rnorm(n_subjects * population$q), n_subjects, population$q)
  b <- sweep(b, 2, population$random_effect_sd, "*")
  noise <- if (error_family == "normal") stats::rnorm(length(rows)) else
    stats::rt(length(rows), df = 3) / sqrt(3)
  y <- as.numeric(population$X[rows, , drop = FALSE] %*% population$beta0) +
    rowSums(population$Z[rows, , drop = FALSE] * b[cluster, , drop = FALSE]) +
    config$errors$marginal_sd * noise
  list(y = y, X = population$X[rows, , drop = FALSE], Z = population$Z[rows, , drop = FALSE],
       cluster_id = cluster, beta0 = population$beta0, source_rows = rows,
       source_subject_types = types, m = sizes, n_clusters = n_subjects,
       random_effects = b, errors = config$errors$marginal_sd * noise,
       error_family = error_family, outcome_source = "declared_matched_simulation")
}

application_population_unit_map_v2 <- function(preprocessing, population_preprocessing, coordinates) {
  pop <- population_preprocessing$gene_preprocessing
  rep <- preprocessing$gene_preprocessing
  M <- diag(length(coordinates))
  dimnames(M) <- list(coordinates, coordinates)
  if (!"(Intercept)" %in% coordinates) stop("Population-unit mapping requires the complete intercept design.")
  modules <- grep("^module__", coordinates, value = TRUE)
  for (coordinate in modules) {
    id <- sub("^module__", "", coordinate)
    if (!id %in% names(pop$scale) || !id %in% names(rep$scale)) stop("Module scale is absent from population mapping.")
    M[coordinate, coordinate] <- pop$scale[id] / rep$scale[id]
    M["(Intercept)", coordinate] <- (pop$center[id] - rep$center[id]) / rep$scale[id]
  }
  if ("centered_time" %in% coordinates) {
    M["(Intercept)", "centered_time"] <- population_preprocessing$clinical$time_center -
      preprocessing$clinical$time_center
  }
  matched <- match(pop$gene_probe_map$gene_symbol, rep$gene_probe_map$gene_symbol)
  if (anyNA(matched)) stop("The reprocessed module gene inventory changed.")
  changes <- sum(pop$gene_probe_map$probe_id != rep$gene_probe_map$probe_id[matched])
  list(matrix = M, changed_gene_probes = changes,
       exact_affine_measurement_map = changes == 0L &&
         identical(preprocessing$clinical$factor_levels, population_preprocessing$clinical$factor_levels) &&
         identical(preprocessing$clinical$numeric_imputation, population_preprocessing$clinical$numeric_imputation),
       measurement_scope = if (changes == 0L) "same_probe_measurements_with_recorded_clinical_preprocessing" else
         "recomputed_probe_measurements_in_population_scale")
}

resample_application_matched_workflow_v2 <- function(population, expression_probe, config,
                                                    error_family, seed) {
  dat <- generate_application_matched_population_data_v2(population, config, error_family,
                                                         config$n_subjects, seed)
  metadata <- population$metadata[dat$source_rows, , drop = FALSE]
  metadata$observed_subject_id <- as.character(metadata$subject_id)
  metadata$observed_sample_id <- rownames(population$metadata)[dat$source_rows]
  metadata$observed_c3 <- metadata$c3
  metadata$subject_id <- sprintf("matched_subject_%06d", dat$cluster_id)
  metadata$sample_id <- paste0(metadata$observed_sample_id, "_draw_", dat$cluster_id)
  rownames(metadata) <- metadata$sample_id
  metadata$c3 <- exp(dat$y)
  metadata$log_c3 <- dat$y
  metadata$outcome_source <- "declared_matched_simulation"
  expression <- expression_probe[, dat$source_rows, drop = FALSE]
  colnames(expression) <- rownames(metadata)
  list(metadata = metadata, expression = expression, generated = dat)
}
