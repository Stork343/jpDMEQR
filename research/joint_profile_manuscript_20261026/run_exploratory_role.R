source("research/score_validity_20261026/load_joint.R")
source("research/joint_profile_manuscript_20261026/fixed_inference.R")
source("research/joint_profile_manuscript_20261026/prediction.R")
args <- commandArgs(trailingOnly = TRUE)
processing <- args[1L]
role <- args[2L]
if (!processing %in% c("reference", "official_R1_R2_normexp_training_quantile", "processed_median_across_mapped_probes")) stop("Unknown exploratory processing route.")
metadata <- readRDS("data/derived/GSE65391/audit_20261022_bootstrap_v1/analysis_cohort.rds")$eligible_metadata
reference_path <- file.path("results/raw/application/sle_complete_20261022_v1_exploratory/preprocessing", paste0(role, ".rds"))
reference <- readRDS(reference_path)
if (processing == "reference") {
  genes <- reference$genes
  clinical <- reference$clinical
  source_path <- reference_path
} else {
  source_path <- file.path("data/derived/GSE65391/processing_20261023_v1", processing,
                           "preprocessing_cache", paste0("exploratory_", role), "prepared_base.rds")
  base <- readRDS(source_path)
  if (!identical(base$sample_ids, rownames(metadata)) || !setequal(base$training_rows, reference$selection_rows)) stop("Exploratory preprocessing roles differ.")
  genes <- base$genes
  clinical <- base$clinical
}
if (!identical(reference$sample_ids, rownames(metadata)) ||
    length(intersect(reference$selection_subjects, reference$inference_subjects))) stop("Invalid exploratory subject roles.")
colnames(genes) <- paste0("gene__", sub("^gene__", "", colnames(genes)))
output <- file.path("results/raw/application_exploratory/joint_profile_20261026", processing, role)
if (dir.exists(output)) stop("This exploratory role already exists.")
dir.create(output, recursive = TRUE)
saveRDS(list(processing = processing, role = role, source_path = source_path,
             selection_subjects = reference$selection_subjects, inference_subjects = reference$inference_subjects,
             additions = 3L, lambda_gamma = 1, started = Sys.time()), file.path(output, "request.rds"))
current <- output
globalCallingHandlers(error = function(condition) saveRDS(condition, file.path(current, "failure.rds")))
jobs <- expand.grid(model = c("A", "B"), tau = c(.25, .5, .75), stringsAsFactors = FALSE)
for (i in seq_len(nrow(jobs))) {
  job <- jobs[i, ]
  current <- file.path(output, sprintf("case_%02d", i))
  dir.create(current)
  saveRDS(job, file.path(current, "request.rds"))
  X <- cbind(clinical[[job$model]]$X, genes)
  selection <- reference$selection_rows
  inference <- reference$inference_rows
  data <- sb_data(metadata$log_c3[selection], X[selection, , drop = FALSE],
                   matrix(1, length(selection), 1L), metadata$subject_id[selection])
  data$penalty_factor <- c(rep(0, ncol(clinical[[job$model]]$X)), rep(1, ncol(genes)))
  prepared <- jp_prepare_prediction_search(data, job$tau, 1)
  path <- pc_orthogonal_path(prepared, job$tau, 0:3, gap_tolerance = 1e-8, max_iterations = 500000L)
  selected <- setdiff(path$selected_path[[4L]], which(data$penalty_factor == 0))
  if (length(selected) != 3L) stop("The exploratory search must select three genes.")
  columns <- c(which(data$penalty_factor == 0), selected)
  selected_names <- colnames(X)[selected]
  fitted <- jp_fit_fixed(metadata$log_c3[inference], X[inference, columns, drop = FALSE],
                          matrix(1, length(inference), 1L), metadata$subject_id[inference], job$tau,
                          coordinates = selected_names)
  table <- fitted$table
  table$gene_id <- sub("^gene__", "", table$coordinate)
  table$model <- job$model
  table$tau <- job$tau
  table$role <- role
  table$processing <- processing
  table$target_scope <- "role_specific_selected_profile"
  saveRDS(list(job = job, table = table, fit = fitted, selected = selected_names,
               selected_path = path$selected_path, search_diagnostics = path$diagnostics,
               source_path = source_path, completed = Sys.time()), file.path(current, "result.rds"))
  write.csv(table, file.path(current, "intervals.csv"), row.names = FALSE)
}
saveRDS(list(completed = nrow(jobs), finished = Sys.time()), file.path(output, "completed.rds"))
