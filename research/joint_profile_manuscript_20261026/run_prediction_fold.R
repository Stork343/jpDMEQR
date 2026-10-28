source("research/score_validity_20261026/load_joint.R")
source_v2_module(".", "application_v2", .GlobalEnv)
source("research/joint_profile_manuscript_20261026/conic_search.R")
source("research/joint_profile_manuscript_20261026/prediction.R")
args <- commandArgs(trailingOnly = TRUE)
processing <- args[1L]
repetition <- as.integer(args[2L])
fold <- as.integer(args[3L])
resume <- length(args) >= 4L && identical(args[4L], "resume")
slope_only <- length(args) >= 4L && identical(args[4L], "slope_only")
if (slope_only && processing != "reference") stop("The slope calculation requires reference processing.")
if (!processing %in% c("reference", "official_R1_R2_normexp_training_quantile", "processed_median_across_mapped_probes") ||
    !repetition %in% 1:5 || !fold %in% 1:5) stop("Invalid prediction role.")
cohort <- readRDS("data/derived/GSE65391/audit_20261022_bootstrap_v1/analysis_cohort.rds")
metadata <- cohort$eligible_metadata
assignments <- cohort$subject_folds
heldout <- assignments$subject_id[assignments$repetition == repetition & assignments$fold == fold]
outer_training <- !metadata$subject_id %in% heldout
inner_roles <- read.csv("data/derived/GSE65391/execution_20261022_v1/inner_subject_folds.csv", stringsAsFactors = FALSE)
inner_roles <- inner_roles[inner_roles$outer_repetition == repetition & inner_roles$outer_fold == fold, ]
if (length(intersect(heldout, inner_roles$subject_id)) ||
    !setequal(inner_roles$subject_id, unique(metadata$subject_id[outer_training]))) stop("Inner and outer roles differ.")
output <- file.path("results/raw/application_prediction/joint_profile_20261026", if (slope_only) "reference_conic_slope" else processing,
                     sprintf("r%02d_f%02d", repetition, fold))
if (dir.exists(output) && !resume) stop("This prediction role already exists.")
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
  saveRDS(list(processing = processing, repetition = repetition, fold = fold,
               selection_subjects = unique(metadata$subject_id[outer_training]),
               heldout_subjects = heldout, inner_roles = inner_roles, started = Sys.time()),
           file.path(output, "request.rds"))
} else {
  prior <- readRDS(file.path(output, "request.rds"))
  if (prior$processing != processing || prior$repetition != repetition || prior$fold != fold ||
      !identical(prior$inner_roles, inner_roles)) stop("The resumed prediction role differs.")
}
current_case <- output
globalCallingHandlers(error = function(condition) saveRDS(list(condition = condition), file.path(current_case, "failure.rds")))
bases <- vector("list", 6L)
if (processing == "reference") {
  expression <- readRDS(cohort$expression_file)[, rownames(metadata), drop = FALSE]
  rule <- yaml::read_yaml("config/module_scoring_20261022.yml")
  modules <- read.csv("config/module_signatures.csv", stringsAsFactors = FALSE)
  members <- read.csv("config/module_score_members_20261022.csv", stringsAsFactors = FALSE)
  annotation <- read.csv(rule$annotation_path, stringsAsFactors = FALSE)
}
for (g in 0:5) {
  train <- if (g == 0L) outer_training else outer_training & !metadata$subject_id %in% inner_roles$subject_id[inner_roles$inner_fold == g]
  if (processing == "reference") {
    design <- prepare_gse65391_prediction_design_v2(metadata, expression, annotation, modules, members,
                                                    rule, train, include_sledai = FALSE, filter_size = 10000L)
    clinical_B <- prepare_clinical_design_split_v2(metadata, train, include_sledai = TRUE)
    gene_columns <- grepl("^gene__", colnames(design$X))
    bases[[g + 1L]] <- list(clinical = list(A = design$clinical, B = clinical_B),
                            module_scores = design$X[, grepl("^module__", colnames(design$X)), drop = FALSE],
                            genes = design$X[, gene_columns, drop = FALSE], training_rows = which(train),
                            sample_ids = rownames(metadata), highdim_preprocessing = design$highdim_preprocessing,
                            gene_preprocessing = design$gene_preprocessing)
  } else {
    suffix <- if (g == 0L) "outer" else sprintf("inner%02d", g)
    path <- file.path("data/derived/GSE65391/processing_20261023_v1", processing,
                       "preprocessing_cache", sprintf("prediction_r%02d_f%02d_%s", repetition, fold, suffix), "prepared_base.rds")
    bases[[g + 1L]] <- readRDS(path)
  }
  if (!identical(bases[[g + 1L]]$sample_ids, rownames(metadata)) ||
      !setequal(bases[[g + 1L]]$training_rows, which(train))) stop("Cached preprocessing subjects differ from the declared role.")
}
if (processing == "reference") rm(expression, design)
specifications <- if (slope_only) {
  data.frame(filter_size = 5000L, random_slope = TRUE)
} else if (processing == "reference") {
  data.frame(filter_size = c(5000L, 2000L, 10000L, 5000L), random_slope = c(FALSE, FALSE, FALSE, TRUE))
} else data.frame(filter_size = 5000L, random_slope = FALSE)
jobs <- expand.grid(model = c("A", "B"), tau = c(.25, .5, .75), specification = seq_len(nrow(specifications)), stringsAsFactors = FALSE)
jobs <- cbind(jobs, specifications[jobs$specification, , drop = FALSE])
write.csv(jobs, file.path(output, "scheduled_cases.csv"), row.names = FALSE)
for (i in seq_len(nrow(jobs))) {
  job <- jobs[i, ]
  current_case <- file.path(output, sprintf("case_%02d", i))
  if (dir.exists(current_case)) {
    if (!resume) stop("Existing prediction case.")
    completed <- list.files(current_case, pattern = "^result[.]rds$", recursive = TRUE)
    if (length(completed) > 1L) stop("Duplicate prediction result.")
    if (length(completed) == 1L) next
    current_case <- file.path(current_case, sprintf("attempt_%02d", 2L + length(list.dirs(current_case, recursive = FALSE))))
  }
  dir.create(current_case, recursive = TRUE)
  saveRDS(job, file.path(current_case, "request.rds"))
  started <- proc.time()[[3L]]
  inner <- lapply(1:5, function(g) {
    train <- bases[[g + 1L]]$training_rows
    valid <- which(metadata$subject_id %in% inner_roles$subject_id[inner_roles$inner_fold == g])
    list(train = jp_prediction_data(metadata, bases[[g + 1L]], train, job$model, job$filter_size, job$random_slope),
         validation = jp_prediction_data(metadata, bases[[g + 1L]], valid, job$model, job$filter_size, job$random_slope))
  })
  outer <- jp_prediction_data(metadata, bases[[1L]], which(outer_training), job$model, job$filter_size, job$random_slope)
  testing <- jp_prediction_data(metadata, bases[[1L]], which(!outer_training), job$model, job$filter_size, job$random_slope)
  fit <- jp_prediction_fit(inner, outer, job$tau,
                           checkpoint = function(record) saveRDS(record, file.path(current_case, "current_stage.rds")))
  prediction <- drop(testing$X %*% fit$beta)
  residual <- testing$y - prediction
  loss <- residual * (job$tau - (residual < 0))
  table <- data.frame(sample_id = rownames(metadata)[!outer_training], subject_id = testing$cluster_id,
                       observed_log_c3 = testing$y, predicted_quantile = prediction, loss = loss,
                       repetition = repetition, fold = fold, model = job$model, tau = job$tau,
                       filter_size = job$filter_size, processing = processing,
                       method = "JOINT-PROFILE", random_slope = job$random_slope)
  fit$runtime <- proc.time()[[3L]] - started
  saveRDS(list(job = job, fit = fit, predictions = table, completed = Sys.time()), file.path(current_case, "result.rds"))
  write.csv(table, file.path(current_case, "predictions.csv"), row.names = FALSE)
  cat("Completed", processing, repetition, fold, i, "of", nrow(jobs), "\n")
}
saveRDS(list(completed = nrow(jobs), finished = Sys.time()), file.path(output, "completed.rds"))
