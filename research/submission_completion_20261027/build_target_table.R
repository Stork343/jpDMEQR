source("research/submission_completion_20261027/load.R")
registry <- submission_registry("main")
records <- lapply(c("B16", "B19"), function(id) {
  config <- config_row_to_list_v2(registry[registry$experiment_id == id, ])
  profile <- readRDS(submission_target_path(config))
  marginal <- readRDS(file.path("results/intermediate/joint_profile_submission_20261027/marginal_targets", paste0(id, ".rds")))
  if (!profile$all_converged || profile$n_population < 100000L || profile$repeats < 4L ||
      !marginal$deterministic_target || marginal$max_kkt_residual > 1e-8 ||
      marginal$quadrature_refinement_difference > 1e-8) stop("Population target checks did not pass.")
  data.frame(experiment_id = id, error = config$error_dist, tau = config$tau,
    generated = make_beta_v2(config$p, config$s, config$signal)[1L],
    profile = profile$beta_star_mc[1L], profile_mcse = profile$target_mc_se_by_coordinate[1L],
    marginal = marginal$beta_marginal[1L], marginal_intercept = marginal$intercept_marginal,
    quadrature_difference = marginal$quadrature_refinement_difference)
})
records <- do.call(rbind, records)
output <- "results/tables/submission_20261028/targets"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
write.csv(records, file.path(output, "population_comparison.csv"), row.names = FALSE)
lines <- c("\\begin{table}[t]",
  "\\caption{Profile and marginal quantile projections under predictor dependent error scales}\\label{tab:conditional-scale-targets}",
  "\\begin{tabularx}{\\textwidth}{@{}Xrrrrrr@{}}", "\\toprule",
  "Error & Quantile & Generated & Profile & MCSE & QR slope & QR intercept \\\\", "\\midrule")
for (i in seq_len(nrow(records))) {
  row <- records[i, ]
  lines <- c(lines, sprintf("%s & %.2f & %.3f & %.3f & %.4f & %.3f & %.3f \\\\",
    if (row$error == "gaussian") "Gaussian" else "$t_3$", row$tau, row$generated,
    row$profile, row$profile_mcse, row$marginal, row$marginal_intercept))
}
lines <- c(lines, "\\bottomrule", "\\end{tabularx}",
  "\\par\\smallskip\\noindent Slope entries concern the first coordinate. Profile targets use four independent populations of 100{,}000 clusters and a half-bandwidth refinement; MCSE refers to that population calculation. Marginal QR includes an unpenalised intercept and uses deterministic quadrature, checked by four independent populations of 100{,}000 clusters. The fitted design has 1000 slope coordinates and eight nonzero generating coefficients.",
  "\\end{table}")
writeLines(lines, "manuscript/generated/tab_conditional_scale_targets.tex")
print(records)
