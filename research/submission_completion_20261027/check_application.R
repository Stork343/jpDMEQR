source("research/score_validity_20261026/load_joint.R")
source_v2_module(".", "application_v2", .GlobalEnv)
source("research/submission_completion_20261027/application.R")
source("research/joint_profile_manuscript_20261026/prediction.R")
cohort <- readRDS("data/derived/GSE65391/audit_20261022_bootstrap_v1/analysis_cohort.rds")
metadata <- cohort$eligible_metadata
heldout <- cohort$subject_folds$subject_id[cohort$subject_folds$repetition == 1L & cohort$subject_folds$fold == 1L]
inner_roles <- read.csv("data/derived/GSE65391/execution_20261022_v1/inner_subject_folds.csv")
inner_roles <- inner_roles[inner_roles$outer_repetition == 1L & inner_roles$outer_fold == 1L, ]
make_data <- function(base, rows, model, tau, transformed) {
  data <- jp_prediction_data(metadata, base, rows, model, 5000L, FALSE)
  if (transformed) {
    clinical <- submission_application_design(base, metadata, model, tau, 0L)$X
    data$X[, seq_len(ncol(clinical))] <- clinical[rows, , drop = FALSE]
  }
  data
}
base <- submission_application_base("reference", 1L, 1L, 0L, metadata)
inner <- lapply(1:5, function(g) {
  local <- submission_application_base("reference", 1L, 1L, g, metadata)
  validation <- which(metadata$subject_id %in% inner_roles$subject_id[inner_roles$inner_fold == g])
  list(train = make_data(local, local$training_rows, "A", .25, TRUE),
       validation = make_data(local, validation, "A", .25, TRUE))
})
outer <- make_data(base, base$training_rows, "A", .25, TRUE)
fitted <- jp_prediction_fit(inner, outer, .25)
stored <- readRDS("results/raw/application_prediction/joint_profile_20261026/reference/r01_f01/case_01/result.rds")$fit
if (!identical(fitted$selected$size, stored$selected$size) || fitted$selected$lambda != stored$selected$lambda ||
    !setequal(names(fitted$beta)[fitted$beta != 0], names(stored$beta)[stored$beta != 0])) stop("Clinical units change the selected joint-profile model.")
original <- make_data(base, base$training_rows, "A", .25, FALSE)
delta <- drop(outer$X %*% fitted$beta - original$X %*% stored$beta)
distance <- sum(outer$weight * delta^2) / outer$n
bound <- 4 * outer$n * (max(fitted$final_gap, 0) + max(stored$final_gap, 0) + 1e-12)
if (distance > bound) stop("Clinical coordinate invariance exceeds the convex certificates.")
checks <- data.frame(check = "joint_profile_clinical_coordinate_invariance", value = distance, bound = bound)
for (tau in c(.25, .5, .75)) {
  design <- submission_application_design(base, metadata, "A", tau, 100L)
  fit <- submission_qgee(metadata, design, seq_len(nrow(metadata)) %in% base$training_rows, tau, 20261027L)
  cat("QGEE quantile", tau, "converged", fit$converged, "penalty", fit$selected_penalty, "\n")
  if (anyNA(fit$beta_prediction)) stop("Unavailable QGEE coefficient.")
  checks <- rbind(checks, data.frame(check = paste0("qgee_returned_tau", tau), value = as.integer(fit$converged), bound = 1))
}
destination <- "results/scratch/submission_application_20261027/checks"
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
write.csv(checks, file.path(destination, "clinical_units.csv"), row.names = FALSE)
saveRDS(list(fitted = fitted, original = stored, distance = distance, bound = bound, checks = checks), file.path(destination, "checked_objects.rds"))
print(checks)
