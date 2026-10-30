source("research/score_validity_20261026/load_joint.R")
source_v2_module(".", "application_v2", .GlobalEnv)
source("research/submission_completion_20261027/dispatch_state.R")
source("research/empirical_recovery_20261030/inputs.R")
input <- empirical_inputs()
base <- empirical_reference_base(input, seq_len(nrow(input$metadata)))
records <- list()
for (model in c("A", "B")) {
  reference <- readRDS(file.path("results/scratch/joint_profile_manuscript_20261026/module_inputs",
    paste0("stages_", model, ".rds")))[[1L]]
  candidate <- cbind(base$clinical[[model]]$X, base$module_scores)
  if (!identical(colnames(candidate), colnames(reference$X))) stop("Module design coordinates differ.")
  difference <- max(abs(candidate - reference$X))
  if (difference > 1e-12) stop("The recovered expression changes the verified module design.")
  records[[model]] <- data.frame(check = paste0("complete_design_", model), difference = difference)
}
roles <- empirical_prediction_roles(input, 1L, 1L)
reference <- empirical_reference_base(input, roles$training[[1L]])
direct <- prepare_gse65391_prediction_design_v2(input$metadata, input$expression,
  input$annotation, input$modules, input$members, input$rule,
  roles$training[[1L]], include_sledai = FALSE, filter_size = 5000L)
genes <- direct$X[, grepl("^gene__", colnames(direct$X)), drop = FALSE]
if (!identical(colnames(reference$genes)[1:5000], colnames(genes)) ||
    max(abs(reference$genes[, 1:5000] - genes)) > 1e-12)
  stop("The shared preprocessing does not reproduce the registered gene filter.")
records[["gene_filter"]] <- data.frame(check = "shared_gene_filter", difference = 0)
output <- "results/preflight/empirical_recovery_20261030"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
write.csv(do.call(rbind, records), file.path(output, "input_checks.csv"), row.names = FALSE)
write.csv(input$fold_audit, file.path(output, "subject_fold_checks.csv"), row.names = FALSE)
write.csv(input$exploratory_audit, file.path(output, "exploratory_role_checks.csv"), row.names = FALSE)
empirical_atomic_json(list(passed = TRUE, subjects = 129L, visits = 867L,
  outer_roles = nrow(input$fold_audit), exploratory_roles = nrow(input$exploratory_audit),
  source = input$source, checked = as.character(Sys.time())),
  file.path(output, "input_manifest.json"))
print(do.call(rbind, records))
