source("research/score_validity_20261026/load_joint.R")
source_v2_module(".", "application_v2", .GlobalEnv)
source("research/submission_completion_20261027/dispatch_state.R")
source("research/empirical_recovery_20261030/inputs.R")
input <- empirical_inputs()
processing <- "processed_median_across_mapped_probes"
median <- empirical_gene_medians(input)
multiple <- names(median$probe_groups)[lengths(median$probe_groups) > 1L]
selected <- unique(c(input$members$member_id[input$members$member_type == "gene_symbol"], head(multiple, 100L)))
errors <- vapply(selected, function(gene) {
  probes <- median$probe_groups[[gene]]
  direct <- apply(input$expression[probes, , drop = FALSE], 2L, stats::median)
  max(abs(direct - median$expression[gene, ]))
}, numeric(1))
if (any(errors > 1e-12)) stop("The median-probe measurements differ from direct aggregation.")
training <- empirical_prediction_roles(input, 1L, 1L)$training[[1L]]
base <- empirical_prepare_base(input, training, processing)
empirical_validate_base(base, input$metadata, training, processing)
mean_error <- max(abs(colMeans(base$module_scores[training, ])))
scale_error <- max(abs(apply(base$module_scores[training, ], 2L, sd) - 1))
if (max(mean_error, scale_error) > 1e-12) stop("Module scaling does not use the declared training observations.")
reference <- empirical_reference_base(input, training)
probe_module_error <- max(abs(base$module_scores[, "module__M5.15"] - reference$module_scores[, "module__M5.15"]))
if (probe_module_error > 1e-12) stop("The probe-defined module measurement differs.")
output <- "results/preflight/empirical_recovery_20261030"
write.csv(data.frame(gene = selected, maximum_difference = errors),
  file.path(output, "median_probe_checks.csv"), row.names = FALSE)
empirical_atomic_json(list(passed = TRUE, processing = processing,
  checked_genes = length(selected), maximum_median_difference = max(errors),
  module_center_error = mean_error, module_scale_error = scale_error,
  probe_module_difference = probe_module_error, source = input$source,
  checked = as.character(Sys.time())), file.path(output, paste0(processing, "_manifest.json")))
cat("Verified median-probe processing for", length(selected), "genes\n")
