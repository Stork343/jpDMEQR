root <- normalizePath(".", mustWork = TRUE)
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", environment())
source_v2_module(root, "application_v2", environment())
args <- commandArgs(trailingOnly = TRUE)
run_id <- if (length(args)) args[1] else "sle_complete_20261022_v1"
version <- if (length(args) > 1L) args[2] else "workload_20261022_v1"
output <- file.path(root, "data", "derived", "GSE65391", version)
if (dir.exists(output)) stop("The full application workload directory already exists.")
config <- yaml::read_yaml(file.path(root, "config", "application_execution_20261022.yml"))
cohort <- readRDS(file.path(root, config$cohort_directory, "analysis_cohort.rds"))
metadata <- cohort$eligible_metadata
if (length(unique(metadata$subject_id)) != 129L || nrow(metadata) != 867L) {
  stop("The complete frozen 129-subject 867-visit cohort is required.")
}
base <- expand.grid(model = config$models, tau = config$quantiles, stringsAsFactors = FALSE)
base$specification <- "reference"
base$h_multiplier <- 1
base$lambda_gamma <- 1
base$random_slope <- FALSE
module <- list(base)
for (value in config$sensitivities$lambda_gamma) {
  item <- base
  item$specification <- paste0("Lambda_", value)
  item$lambda_gamma <- value
  module[[length(module) + 1L]] <- item
}
for (value in config$sensitivities$h_multipliers) {
  item <- base
  item$specification <- paste0("h_", value)
  item$h_multiplier <- value
  module[[length(module) + 1L]] <- item
}
item <- base
item$specification <- "random_slope"
item$random_slope <- TRUE
module[[length(module) + 1L]] <- item
curvature <- cohort$audit$curvature$diagnostics
eligible <- curvature$h_multiplier[curvature$tau == .1 & curvature$pass]
if (!setequal(eligible, c(1, 1.25))) stop("The frozen tau .10 curvature eligibility changed.")
item <- expand.grid(model = config$models, tau = .1, h_multiplier = eligible, stringsAsFactors = FALSE)
item$specification <- paste0("tau10_h_", item$h_multiplier)
item$lambda_gamma <- 1
item$random_slope <- FALSE
module[[length(module) + 1L]] <- item[names(base)]
module <- do.call(rbind, module)
module$job_id <- seq_len(nrow(module))
module$run_id <- sprintf("%s_module_%03d_%s_tau%02d", run_id, module$job_id, module$model, 100 * module$tau)
module$expected_full_fits <- 1L
module$expected_leave_one_fits <- 129L
module$module_count <- 5L
module$coordinate_count <- ifelse(module$model == "A", 16L, 17L)
module$entrypoint <- "scripts/application/03_fit_lowdim_GSE65391.R"
module$expected_manifest <- file.path("results", "raw", "application", module$run_id, "completed_manifest.json")
module_commands <- lapply(seq_len(nrow(module)), function(i) {
  job <- module[i, ]
  list(job = job, script = job$entrypoint, arguments = c(
    paste0("--cohort-dir=", file.path(root, config$cohort_directory)), paste0("--models=", job$model),
    paste0("--quantiles=", job$tau), paste0("--h-multiplier=", job$h_multiplier),
    paste0("--lambda-gamma=", job$lambda_gamma),
    paste0("--random-slope=", tolower(as.character(job$random_slope))), paste0("--run-id=", job$run_id)
  ))
})
specifications <- data.frame(
  method = c("PROFILE-INITIAL", "POOLED-QR-LASSO", "QGEE-SCAD", "LQMM",
             "PROFILE-INITIAL", "PROFILE-INITIAL", "PROFILE-INITIAL"),
  filter_size = c(5000L, 5000L, 5000L, 5000L, 2000L, 10000L, 5000L),
  random_slope = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
  stringsAsFactors = FALSE
)
validation <- expand.grid(repetition = seq_len(5L), fold = seq_len(5L), model = config$models,
                           tau = config$quantiles, specification_id = seq_len(nrow(specifications)),
                           stringsAsFactors = FALSE)
validation <- cbind(validation, specifications[validation$specification_id, , drop = FALSE])
rownames(validation) <- NULL
validation$job_id <- seq_len(nrow(validation))
validation$run_id <- sprintf("%s_cv_%04d_%s_%s_r%d_f%d_tau%02d_filter%d_q%d", run_id,
                              validation$job_id, validation$method, validation$model,
                              validation$repetition, validation$fold, 100 * validation$tau,
                              validation$filter_size, 1L + as.integer(validation$random_slope))
candidates <- c("PROFILE-INITIAL" = 108L, "POOLED-QR-LASSO" = 12L, "QGEE-SCAD" = 3L, "LQMM" = 1L)
validation$candidates <- unname(candidates[validation$method])
validation$inner_folds <- 5L
validation$expected_inner_fits <- validation$candidates * validation$inner_folds
validation$expected_outer_fits <- 1L
validation$entrypoint <- "scripts/application/12_validate_GSE65391.R"
validation$expected_manifest <- file.path("results", "raw", "application_validation", validation$run_id, "completed_manifest.json")
validation_commands <- lapply(seq_len(nrow(validation)), function(i) {
  job <- validation[i, ]
  list(job = job, script = job$entrypoint, arguments = c(
    paste0("--repetition=", job$repetition), paste0("--fold=", job$fold), paste0("--model=", job$model),
    paste0("--tau=", job$tau), paste0("--features=", job$filter_size), paste0("--method=", job$method),
    paste0("--random-slope=", tolower(as.character(job$random_slope))), paste0("--run-id=", job$run_id)
  ))
})
roles <- read.csv(file.path(root, "data", "derived", "GSE65391", "execution_20261022_v1",
                            "exploratory_subject_roles.csv"), stringsAsFactors = FALSE)
audit <- validate_exploratory_subject_roles_v2(metadata, roles)
exploratory <- expand.grid(model = config$models, tau = config$quantiles,
                            role_id = audit$role_id, stringsAsFactors = FALSE)
exploratory <- merge(exploratory, audit[c("role_id", "split", "swapped")], by = "role_id", sort = FALSE)
exploratory <- exploratory[order(exploratory$split, exploratory$swapped, exploratory$model, exploratory$tau), ]
rownames(exploratory) <- NULL
exploratory$job_id <- seq_len(nrow(exploratory))
exploratory$run_id <- paste0(run_id, "_exploratory")
exploratory$filter_size <- 5000L
exploratory$selected_genes <- 3L
exploratory$clinical_count <- ifelse(exploratory$model == "A", 11L, 12L)
exploratory$entrypoint <- "scripts/application/09_fit_exploratory_pair_GSE65391.R"
exploratory$expected_manifest <- file.path("results", "raw", "application", exploratory$run_id,
                                           sprintf("job%03d", exploratory$job_id), "completed_manifest.json")
exploratory_commands <- lapply(seq_len(nrow(exploratory)), function(i) {
  job <- exploratory[i, ]
  list(job = job, script = job$entrypoint, arguments = c(
    paste0("--job-index=", job$job_id), paste0("--run-id=", job$run_id)
  ))
})
if (nrow(module) != 40L || nrow(validation) != 1050L || nrow(exploratory) != 120L ||
    anyDuplicated(module$run_id) || anyDuplicated(validation$run_id)) {
  stop("The complete predeclared application workload matrix changed.")
}
dir.create(output, recursive = TRUE)
write.csv(module, file.path(output, "module_jobs.csv"), row.names = FALSE)
write.csv(validation, file.path(output, "validation_jobs.csv"), row.names = FALSE)
write.csv(exploratory, file.path(output, "exploratory_jobs.csv"), row.names = FALSE)
saveRDS(list(module = module_commands, validation = validation_commands, exploratory = exploratory_commands),
        file.path(output, "commands.rds"))
file.copy(file.path(root, "config", "application_execution_20261022.yml"), file.path(output, "application_execution.yml"))
jsonlite::write_json(list(
  run_id = run_id, n_subjects = 129L, n_visits = 867L, module_count = 5L,
  module_jobs = nrow(module), full_workflow_fits = sum(module$expected_full_fits + module$expected_leave_one_fits),
  validation_outer_jobs = nrow(validation), validation_inner_fits = sum(validation$expected_inner_fits),
  exploratory_jobs = nrow(exploratory), exploratory_role_count = nrow(audit),
  all_job_count = nrow(module) + nrow(validation) + nrow(exploratory),
  raw_route = "official_R1_R2_preprocessing_to_follow_quantitative_processed_route",
  required_gate = "actual_matching_paired_application_method_checks",
  construction_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  created_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
), file.path(output, "workload_manifest.json"), auto_unbox = TRUE, pretty = TRUE)
