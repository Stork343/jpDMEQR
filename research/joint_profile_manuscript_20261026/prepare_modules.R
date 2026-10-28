source("research/score_validity_20261026/load_joint.R")
source_v2_module(".", "application_v2", .GlobalEnv)
root <- "results/scratch/joint_profile_manuscript_20261026/module_inputs"
if (dir.exists(root)) stop("The module input directory already exists.")
cohort_path <- "data/derived/GSE65391/audit_20261022_bootstrap_v1/analysis_cohort.rds"
cohort <- readRDS(cohort_path)
modules <- validate_gse65391_fit_inputs_v2(cohort$audit, "config/module_signatures.csv")
members <- read_frozen_module_score_members_v2("config/module_score_members_20261022.csv")
rule <- yaml::read_yaml("config/module_scoring_20261022.yml")
annotation <- read.csv(rule$annotation_path, stringsAsFactors = FALSE)
metadata <- cohort$eligible_metadata
expression_path <- "data/derived/GSE65391/module_measurements_20261022_v1/observed_module_probe_expression.rds"
expression <- readRDS(expression_path)
if (!identical(colnames(expression), rownames(metadata)) || length(unique(metadata$subject_id)) != 129L || nrow(metadata) != 867L) stop("The fixed cohort or measurement order differs.")
checks <- read.csv("results/scratch/joint_profile_manuscript_20261026/fixed_checks.csv")
if (nrow(checks) != 18L || any(checks$estimate_error > 2e-6) || any(checks$se_error > 2e-6)) stop("The fixed inference checks are incomplete.")
dir.create(root, recursive = TRUE)
subjects <- unique(as.character(metadata$subject_id))
saveRDS(list(cohort_path = cohort_path, expression_path = expression_path, module_rule = rule,
             module_members = members, subjects = subjects, sample_ids = rownames(metadata),
             created = Sys.time()), file.path(root, "request.rds"))
for (model in c("A", "B")) {
  stages <- vector("list", length(subjects) + 1L)
  expected <- NULL
  for (i in 0:length(subjects)) {
    keep <- if (i == 0L) rep(TRUE, nrow(metadata)) else metadata$subject_id != subjects[i]
    stage <- prepare_gse65391_profile_workflow_stage_v2(metadata[keep, , drop = FALSE],
                                                       expression[, keep, drop = FALSE], annotation,
                                                       modules, members, rule, include_sledai = model == "B",
                                                       expected_coordinates = expected)
    if (i == 0L) expected <- colnames(stage$X)
    stage$time <- metadata$cumulative_time[keep] - stage$preprocessing$clinical$time_center
    stage$omitted_subject <- if (i == 0L) NA_character_ else subjects[i]
    stage$preprocessing <- compact_application_preprocessing_v2(stage$preprocessing)
    stages[[i + 1L]] <- stage
  }
  saveRDS(stages, file.path(root, paste0("stages_", model, ".rds")))
  cat("Prepared", model, "with", length(stages), "workflow designs\n")
}
base <- expand.grid(model = c("A", "B"), tau = c(.25, .5, .75), stringsAsFactors = FALSE)
base$specification <- "reference"
base$h_multiplier <- 1
base$lambda_gamma <- 1
base$random_slope <- FALSE
jobs <- list(base)
for (lambda in c(.5, 2)) {
  setting <- base
  setting$lambda_gamma <- lambda
  setting$specification <- paste0("Lambda_", lambda)
  jobs[[length(jobs) + 1L]] <- setting
}
for (multiplier in c(.75, 1.25)) {
  setting <- base
  setting$h_multiplier <- multiplier
  setting$specification <- paste0("h_", multiplier)
  jobs[[length(jobs) + 1L]] <- setting
}
setting <- base
setting$random_slope <- TRUE
setting$specification <- "random_slope"
jobs[[length(jobs) + 1L]] <- setting
eligible <- cohort$audit$curvature$diagnostics
eligible <- eligible$h_multiplier[eligible$tau == .1 & eligible$pass]
if (!setequal(eligible, c(1, 1.25))) stop("The lower quantile eligibility changed.")
setting <- expand.grid(model = c("A", "B"), tau = .1, h_multiplier = eligible, stringsAsFactors = FALSE)
setting$specification <- paste0("tau10_h_", setting$h_multiplier)
setting$lambda_gamma <- 1
setting$random_slope <- FALSE
jobs[[length(jobs) + 1L]] <- setting[names(base)]
jobs <- do.call(rbind, jobs)
jobs$job_id <- seq_len(nrow(jobs))
stopifnot(nrow(jobs) == 40L)
write.csv(jobs, file.path(root, "jobs.csv"), row.names = FALSE)
