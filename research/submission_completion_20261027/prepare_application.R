source("research/score_validity_20261026/load_joint.R")
source_v2_module(".", "application_v2", .GlobalEnv)
args <- commandArgs(trailingOnly = TRUE)
repetition <- as.integer(args[1L])
fold <- as.integer(args[2L])
if (!repetition %in% 1:5 || !fold %in% 1:5) stop("Invalid outer role.")
cohort <- readRDS("data/derived/GSE65391/audit_20261022_bootstrap_v1/analysis_cohort.rds")
metadata <- cohort$eligible_metadata
heldout <- cohort$subject_folds$subject_id[cohort$subject_folds$repetition == repetition & cohort$subject_folds$fold == fold]
outer_training <- !metadata$subject_id %in% heldout
inner <- read.csv("data/derived/GSE65391/execution_20261022_v1/inner_subject_folds.csv")
inner <- inner[inner$outer_repetition == repetition & inner$outer_fold == fold, ]
expression <- readRDS(cohort$expression_file)[, rownames(metadata), drop = FALSE]
rule <- yaml::read_yaml("config/module_scoring_20261022.yml")
modules <- read.csv("config/module_signatures.csv")
members <- read.csv("config/module_score_members_20261022.csv")
annotation <- read.csv(rule$annotation_path)
for (g in 0:5) {
  train <- if (g == 0L) outer_training else outer_training & !metadata$subject_id %in% inner$subject_id[inner$inner_fold == g]
  suffix <- if (g == 0L) "outer" else sprintf("inner%02d", g)
  destination <- file.path("results/scratch/submission_application_20261027/preprocessing", sprintf("prediction_r%02d_f%02d_%s", repetition, fold, suffix))
  if (dir.exists(destination)) stop("The reference preprocessing role already exists.")
  dir.create(destination, recursive = TRUE)
  design <- prepare_gse65391_prediction_design_v2(metadata, expression, annotation, modules, members, rule, train,
                                                 include_sledai = FALSE, filter_size = 5000L)
  base <- list(clinical = list(A = design$clinical, B = prepare_clinical_design_split_v2(metadata, train, include_sledai = TRUE)),
               module_scores = design$X[, grepl("^module__", colnames(design$X)), drop = FALSE],
               genes = design$X[, grepl("^gene__", colnames(design$X)), drop = FALSE],
               training_rows = which(train), training_subjects = unique(metadata$subject_id[train]),
               sample_ids = rownames(metadata), highdim_preprocessing = design$highdim_preprocessing,
               gene_preprocessing = design$gene_preprocessing)
  saveRDS(base, file.path(destination, "prepared_base.rds"), compress = FALSE)
}
cat("Prepared reference role", repetition, fold, "\n")
