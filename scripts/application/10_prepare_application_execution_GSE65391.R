root <- normalizePath(".", mustWork = TRUE)
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", envir = environment())
source_v2_module(root, "application_v2", envir = environment())
config <- yaml::read_yaml(file.path(root, "config", "application_execution_20261022.yml"))
matched <- yaml::read_yaml(file.path(root, "config", "application_matched_prospective_20261022.yml"))
output <- file.path(root, "data", "derived", "GSE65391", "execution_20261022_v1")
if (dir.exists(output)) stop("The prospective application execution definition already exists.")
obj <- readRDS(file.path(root, config$subject_design_directory, "analysis_inputs.rds"))
metadata <- obj$eligible_metadata
outer <- obj$subject_folds
validate_gse65391_subject_folds_v2(metadata, outer)
inner <- list()
for (r in seq_len(5L)) for (f in seq_len(5L)) {
  heldout <- outer$subject_id[outer$repetition == r & outer$fold == f]
  training <- metadata[!metadata$subject_id %in% heldout, , drop = FALSE]
  fold <- make_repeated_subject_folds_v2(
    training, rep(TRUE, nrow(training)), k = config$inner_folds, repetitions = 1L,
    seed = config$base_seed + r * 100000L + f * 1009L
  )
  fold$outer_repetition <- r
  fold$outer_fold <- f
  names(fold)[names(fold) == "fold"] <- "inner_fold"
  fold$repetition <- NULL
  if (length(intersect(fold$subject_id, heldout))) stop("Outer held-out subjects entered an inner fold.")
  inner[[length(inner) + 1L]] <- fold
}
inner <- do.call(rbind, inner)
subjects <- sort(unique(as.character(metadata$subject_id)))
exploratory <- list()
for (split in seq_len(config$exploratory$subject_splits)) {
  set.seed(config$base_seed + split * 7919L)
  selection <- sample(subjects, floor(length(subjects) * config$exploratory$selection_fraction))
  for (swap in c(FALSE, TRUE)) {
    group <- if (swap) setdiff(subjects, selection) else selection
    exploratory[[length(exploratory) + 1L]] <- data.frame(
      split = split, swapped = swap, subject_id = subjects,
      role = ifelse(subjects %in% group, "selection", "inference"), stringsAsFactors = FALSE
    )
  }
}
exploratory <- do.call(rbind, exploratory)
dir.create(output, recursive = TRUE)
write.csv(inner, file.path(output, "inner_subject_folds.csv"), row.names = FALSE)
write.csv(exploratory, file.path(output, "exploratory_subject_roles.csv"), row.names = FALSE)
rule <- yaml::read_yaml(file.path(root, "config", "module_scoring_20261022.yml"))
members <- read_frozen_module_score_members_v2(file.path(root, rule$score_membership_file))
signatures <- read_frozen_module_definitions_v2(file.path(root, rule$membership_file))
annotation <- read.csv(file.path(root, rule$annotation_path), stringsAsFactors = FALSE)
expression <- readRDS(obj$expression_file)[, rownames(metadata), drop = FALSE]
for (model in config$models) for (q in matched$random_effect_dimensions) {
  population <- prepare_application_matched_population_v2(
    metadata, expression, annotation, signatures, members, rule, matched, model, q
  )
  identity_map <- application_population_unit_map_v2(
    population$preprocessing, population$preprocessing, colnames(population$X)
  )
  if (max(abs(identity_map$matrix - diag(ncol(population$X)))) > 1e-12) {
    stop("The same population preprocessing must produce the identity unit map.")
  }
  saveRDS(population, file.path(output, sprintf("matched_population_model_%s_q%d.rds", model, q)))
}
registry <- expand.grid(model = matched$models, q = matched$random_effect_dimensions,
                         error_family = matched$error_families, tau = matched$quantiles,
                         stringsAsFactors = FALSE)
registry$case_id <- sprintf("AM_%s_q%d_%s_tau%02d", registry$model, registry$q,
                             registry$error_family, round(100 * registry$tau))
registry$replications <- matched$replications
write.csv(registry, file.path(output, "matched_registry.csv"), row.names = FALSE)
file.copy(file.path(root, "config", "application_execution_20261022.yml"), output)
file.copy(file.path(root, "config", "application_matched_prospective_20261022.yml"), output)
jsonlite::write_json(list(
  n_subjects = length(subjects), n_visits = nrow(metadata), outer_folds = 25L,
  inner_assignment_rows = nrow(inner), exploratory_assignment_rows = nrow(exploratory),
  matched_cases = nrow(registry), matched_replications_per_case = matched$replications,
  inner_md5 = unname(tools::md5sum(file.path(output, "inner_subject_folds.csv"))),
  exploratory_md5 = unname(tools::md5sum(file.path(output, "exploratory_subject_roles.csv"))),
  complete_association_status = "awaiting_matching_method_gate",
  created_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
), file.path(output, "execution_manifest.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
print(registry)
