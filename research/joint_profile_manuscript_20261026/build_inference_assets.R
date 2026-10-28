library(ggplot2)
library(dplyr)

output <- "results/tables/joint_profile_manuscript_20261026"
generated <- "manuscript/generated/joint_profile"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
dir.create(generated, recursive = TRUE, showWarnings = FALSE)
fmt <- function(x, digits = 3L) formatC(ifelse(abs(x) < .5 * 10^(-digits), 0, x), digits = digits, format = "f")
interval <- function(lo, hi) paste0("[", fmt(lo), ", ", fmt(hi), "]")
tex_table <- function(rows, header, caption, label, file, long = FALSE, note = NULL) {
  alignment <- paste0("l", paste(rep("r", ncol(rows) - 1L), collapse = ""))
  body <- apply(rows, 1L, function(row) paste0(paste(row, collapse = " & "), " \\\\"))
  if (long) {
    lines <- c(paste0("\\begin{longtable}{", alignment, "}"),
               paste0("\\caption{", caption, "}\\label{", label, "}\\\\"), "\\toprule",
               paste0(paste(header, collapse = " & "), " \\\\"), "\\midrule\\endfirsthead",
               "\\toprule", paste0(paste(header, collapse = " & "), " \\\\"), "\\midrule\\endhead",
               body, "\\bottomrule", "\\end{longtable}")
  } else {
    lines <- c("\\begin{table}[t]", paste0("\\caption{", caption, "}\\label{", label, "}"),
               "\\centering\\small", paste0("\\begin{tabular}{", alignment, "}"), "\\toprule",
               paste0(paste(header, collapse = " & "), " \\\\"), "\\midrule", body,
               "\\bottomrule", "\\end{tabular}", if (!is.null(note)) paste0("\\par\\smallskip\\parbox{\\textwidth}{\\footnotesize ", note, "}"), "\\end{table}")
  }
  writeLines(lines, file.path(generated, paste0(file, ".tex")))
}
macros <- character()
macro <- function(name, value) {
  if (name %in% names(macros)) stop("Duplicate manuscript value.")
  macros[name] <<- paste0("\\newcommand{\\", name, "}{", value, "}")
}

registry <- read.csv("research/score_validity_20261026/evaluation.csv")
full <- read.csv("results/tables/score_validity_20261026/design_independent/coordinates.csv")
ablation_files <- list.files("results/scratch/joint_profile_manuscript_20261026/ablation", "^result[.]rds$", recursive = TRUE, full.names = TRUE)
if (length(ablation_files) != sum(registry$replications)) stop("Incomplete paired simulation ablation.")
ablated <- bind_rows(lapply(ablation_files, function(path) readRDS(path)$table))
if (nrow(full) != 4L * sum(registry$replications) || nrow(ablated) != nrow(full)) stop("Incomplete coordinate results.")
paired <- full %>% inner_join(ablated, by = c("case_id", "replicate", "coordinate", "truth"), suffix = c("", "_ablated"))
if (nrow(paired) != nrow(full) || anyDuplicated(paired[c("case_id", "replicate", "coordinate")])) stop("Paired simulation keys differ.")
paired <- paired %>% mutate(length = upper - lower, length_ablated = upper_ablated - lower_ablated,
                            difference = as.numeric(covered) - as.numeric(covered_ablated)) %>%
  left_join(registry, by = "case_id")
paired$group <- with(paired, ifelse(grepl("^WEAK", case_id), "Weak signals",
  ifelse(grepl("^T3", case_id), "$t_3$ errors", ifelse(rho_x == .8, "Correlation 0.8",
  ifelse(p > 500, "Larger dimension", ifelse(n_clusters == 200, "Gaussian n200", "Gaussian n400"))))))
replicate_means <- paired %>% group_by(case_id, group, type, replicate) %>%
  summarise(coverage = mean(covered), coverage_ablated = mean(covered_ablated),
            difference = mean(difference), length = mean(length), length_ablated = mean(length_ablated), .groups = "drop")
scenario <- replicate_means %>% group_by(case_id, group, type) %>%
  summarise(replications = n(), mcse = sd(coverage) / sqrt(n()), coverage = mean(coverage),
            ablation_mcse = sd(coverage_ablated) / sqrt(n()), coverage_ablated = mean(coverage_ablated),
            difference_mcse = sd(difference) / sqrt(n()), difference = mean(difference),
            length = mean(length), length_ablated = mean(length_ablated), .groups = "drop")
grouped <- scenario %>% group_by(group, type) %>% summarise(settings = n(),
  datasets = sum(replications), coverage = mean(coverage), mcse = sqrt(sum(mcse^2)) / n(),
  coverage_ablated = mean(coverage_ablated), ablation_mcse = sqrt(sum(ablation_mcse^2)) / n(),
  difference = mean(difference), difference_mcse = sqrt(sum(difference_mcse^2)) / n(),
  length = mean(length), length_ablated = mean(length_ablated), .groups = "drop")
overall <- scenario %>% group_by(type) %>% summarise(coverage = mean(coverage),
  mcse = sqrt(sum(mcse^2)) / n(), coverage_ablated = mean(coverage_ablated),
  difference = mean(difference), difference_mcse = sqrt(sum(difference_mcse^2)) / n(),
  length = mean(length), length_ablated = mean(length_ablated), .groups = "drop")
coordinate <- paired %>% group_by(case_id, coordinate, type) %>% summarise(
  replications = n(), coverage = mean(covered), covered = sum(covered),
  bias = mean(beta_tilde - truth), bias_mcse = sd(beta_tilde - truth) / sqrt(n()),
  sampling_sd = sd(beta_tilde), mean_se = mean(se), interval_length = mean(length), .groups = "drop")
binomial <- t(vapply(seq_len(nrow(coordinate)), function(i) binom.test(coordinate$covered[i], coordinate$replications[i])$conf.int, numeric(2)))
coordinate$binomial_lower <- binomial[, 1L]
coordinate$binomial_upper <- binomial[, 2L]
coordinate$coverage_mcse <- sqrt(coordinate$coverage * (1 - coordinate$coverage) / coordinate$replications)
if (any(!is.finite(scenario$mcse)) || any(!is.finite(coordinate$coverage_mcse)) ||
    any(coordinate$coverage < 0 | coordinate$coverage > 1)) stop("Invalid Monte Carlo summaries.")
write.csv(paired, file.path(output, "paired_simulation_coordinates.csv"), row.names = FALSE)
write.csv(scenario, file.path(output, "simulation_scenarios.csv"), row.names = FALSE)
write.csv(grouped, file.path(output, "simulation_groups.csv"), row.names = FALSE)
write.csv(overall, file.path(output, "simulation_overall.csv"), row.names = FALSE)
write.csv(coordinate, file.path(output, "simulation_coordinate_intervals.csv"), row.names = FALSE)
format_sim <- function(x) data.frame(Setting = ifelse(x$group == "$t_3$ errors", x$group, gsub("_", "\\\\_", x$group)),
  Coordinate = ifelse(x$type == "nonzero", "Active", "Null"),
  Full = paste0(fmt(x$coverage), " (", fmt(x$mcse), ")"),
  Ablated = paste0(fmt(x$coverage_ablated), " (", fmt(x$ablation_mcse), ")"),
  Difference = paste0(fmt(x$difference), " (", fmt(x$difference_mcse), ")"),
  Length = paste0(fmt(x$length), " / ", fmt(x$length_ablated)))
tex_table(format_sim(grouped), c("Setting", "Coordinate", "Full coverage", "Ablation", "Difference", "Length"),
  "Coverage and interval length with paired removal of design adjustment", "tab:joint-simulation",
  "simulation_groups", note = "Parentheses give Monte Carlo standard errors. Differences use the same generated datasets. Length lists the full procedure followed by the ablation. Each dataset contributes two active and two null coordinates; standard errors use averages within datasets and retain the scenario strata.")
scenario$group <- scenario$case_id
tex_table(format_sim(scenario), c("Setting", "Coordinate", "Full coverage", "Ablation", "Difference", "Length"),
  "Scenario results for the paired simulation comparison", "tab:joint-simulation-all", "simulation_scenarios", long = TRUE)
tex_table(data.frame(Setting = gsub("_", "\\\\_", registry$case_id), n = registry$n_clusters, p = registry$p,
                    q = registry$q, tau = fmt(registry$tau, 2), Error = registry$error_dist,
                    Signal = fmt(registry$signal, 2), Correlation = fmt(registry$rho_x, 1)),
          c("Setting", "$n$", "$p$", "$q$", "$\\tau$", "Error", "Signal", "$\\rho_X$"),
          "Simulation configurations with sixty independent datasets per row", "tab:joint-registry", "simulation_registry", long = TRUE)
tex_table(data.frame(Setting = gsub("_", "\\\\_", coordinate$case_id), Coordinate = coordinate$coordinate,
                    Coverage = fmt(coordinate$coverage), MCSE = fmt(coordinate$coverage_mcse),
                    Interval = interval(coordinate$binomial_lower, coordinate$binomial_upper),
                    Ratio = fmt(coordinate$mean_se / coordinate$sampling_sd)),
          c("Setting", "Coordinate", "Coverage", "MCSE", "Binomial interval", "SE/SD"),
          "Coordinate coverage with exact binomial intervals and standard error ratios", "tab:joint-coordinate", "simulation_coordinates", long = TRUE)
for (type in c("nonzero", "zero")) {
  row <- overall[overall$type == type, ]
  prefix <- if (type == "nonzero") "JointActive" else "JointNull"
  macro(paste0(prefix, "Coverage"), fmt(row$coverage))
  macro(paste0(prefix, "MCSE"), fmt(row$mcse))
  macro(paste0(prefix, "AblationCoverage"), fmt(row$coverage_ablated))
  macro(paste0(prefix, "CoverageDifference"), fmt(row$difference))
  macro(paste0(prefix, "DifferenceMCSE"), fmt(row$difference_mcse))
  macro(paste0(prefix, "Length"), fmt(row$length))
  macro(paste0(prefix, "AblationLength"), fmt(row$length_ablated))
}
macro("JointDatasets", nrow(full) / 4L)
macro("JointIntervals", nrow(full))
macro("JointAblationShiftCount", sum(abs(paired$target_shift) > 1e-8))
macro("JointAblationMaxShift", fmt(max(abs(paired$target_shift))))
weak <- grouped[grouped$group == "Weak signals" & grouped$type == "nonzero", ]
macro("JointWeakCoverage", fmt(weak$coverage))
macro("JointWeakAblationCoverage", fmt(weak$coverage_ablated))
macro("JointWeakDifference", fmt(weak$difference))
macro("JointWeakDifferenceMCSE", fmt(weak$difference_mcse))
target_audit <- read.csv("results/tables/score_validity_20261026/design_independent/selected_target_audit.csv")
if (nrow(target_audit) != 2L * nrow(full) || any(abs(target_audit$target_shift) > 1e-8)) stop("The full method target comparison differs from the reported identity.")
write.csv(target_audit, file.path(output, "selected_target_audit.csv"), row.names = FALSE)
scaling <- read.csv("results/tables/score_validity_20261026/design_independent/scaling.csv")
tex_table(data.frame(p = scaling$p, Repetitions = scaling$replications, Median = fmt(scaling$median),
                    P90 = fmt(scaling$p90), Maximum = fmt(scaling$maximum),
                    Inference = fmt(scaling$median_inference)),
          c("$p$", "Repetitions", "Median", "90th percentile", "Maximum", "Inference"),
          "Seconds for four requested coordinates at two hundred clusters", "tab:joint-runtime", "runtime")
macro("JointRuntimeFiveHundred", fmt(scaling$median[scaling$p == 500]))
macro("JointRuntimeTenThousand", fmt(scaling$median[scaling$p == 10000]))
orders <- read.csv(file.path(output, "population_orders.csv"))
macro("JointSingleMinimumOrder", fmt(min(orders$single_order)))
macro("JointSingleMaximumOrder", fmt(max(orders$single_order)))
macro("JointPairMinimumOrder", fmt(min(orders$paired_order)))
macro("JointPairMaximumOrder", fmt(max(orders$paired_order)))
write.csv(scaling, file.path(output, "runtime.csv"), row.names = FALSE)

module_files <- c(list.files("results/raw/application/joint_profile_20261026", "^result[.]rds$", recursive = TRUE, full.names = TRUE),
                  list.files("results/raw/application/joint_profile_processing_20261026", "^result[.]rds$", recursive = TRUE, full.names = TRUE))
if (length(module_files) != 52L) stop("Incomplete module workflows.")
macro("JointModuleWorkflows", length(module_files))
macro("JointModulePairedFits", length(module_files) * 130L)
modules <- bind_rows(lapply(module_files, function(path) {
  result <- readRDS(path)
  if (nrow(result$diagnostics) != 130L) stop("Incomplete subject deletion workflow.")
  table <- result$table[grepl("^module__", result$table$coordinate), ]
  if (nrow(table) != 5L) stop("Incomplete module family.")
  if (!"processing" %in% names(table)) table$processing <- "reference"
  if (!"specification" %in% names(table)) {
    table$specification <- "reference"
    table$lambda_gamma <- 1
    table$h_multiplier <- 1
    table$random_slope <- FALSE
  }
  table$source_file <- path
  table
}))
modules$family_size <- ifelse(modules$tau == .1, 10L, 30L)
critical <- qnorm(1 - .05 / (2 * modules$family_size))
modules$simultaneous_lower <- modules$beta_hat - critical * modules$se
modules$simultaneous_upper <- modules$beta_hat + critical * modules$se
modules$workflow_simultaneous_lower <- modules$beta_hat - critical * modules$workflow_se
modules$workflow_simultaneous_upper <- modules$beta_hat + critical * modules$workflow_se
modules$module <- sub("^module__", "", modules$coordinate)
modules$score_excludes_zero <- modules$simultaneous_lower * modules$simultaneous_upper > 0
modules$workflow_excludes_zero <- modules$workflow_simultaneous_lower * modules$workflow_simultaneous_upper > 0
if (any(!is.finite(as.matrix(modules[c("beta_hat", "se", "workflow_se", "simultaneous_lower", "simultaneous_upper")])))) stop("Unavailable module interval values.")
write.csv(modules, file.path(output, "module_intervals.csv"), row.names = FALSE)
reference <- modules %>% filter(processing == "reference", specification == "reference") %>% arrange(model, tau, module)
if (nrow(reference) != 30L) stop("Reference module family differs.")
tex_table(data.frame(Module = reference$module, Model = reference$model, Quantile = fmt(reference$tau, 2),
                    Estimate = fmt(reference$beta_hat), Score = interval(reference$simultaneous_lower, reference$simultaneous_upper),
                    Workflow = interval(reference$workflow_simultaneous_lower, reference$workflow_simultaneous_upper)),
          c("Module", "Model", "$\\tau$", "Estimate", "Score interval", "Workflow interval"),
          "The complete reference module family", "tab:joint-modules", "module_intervals", long = TRUE)
processing_modules <- modules %>% filter(processing != "reference") %>% arrange(processing, model, tau, module)
processing_modules$route <- ifelse(grepl("official", processing_modules$processing), "Raw arrays", "Median probes")
tex_table(data.frame(Route = processing_modules$route, Module = processing_modules$module,
                    Model = processing_modules$model, Quantile = fmt(processing_modules$tau, 2),
                    Estimate = fmt(processing_modules$beta_hat), Score = interval(processing_modules$simultaneous_lower, processing_modules$simultaneous_upper),
                    Workflow = interval(processing_modules$workflow_simultaneous_lower, processing_modules$workflow_simultaneous_upper)),
          c("Processing", "Module", "Model", "$\\tau$", "Estimate", "Score interval", "Workflow interval"),
          "Module intervals under alternative expression processing", "tab:joint-processing-modules", "processing_modules", long = TRUE)
sensitivity <- modules %>% filter(processing == "reference", specification != "reference") %>%
  group_by(specification, lambda_gamma, h_multiplier, random_slope) %>% summarise(
    workflows = n_distinct(source_file), intervals = n(), score_exclusions = sum(score_excludes_zero),
    workflow_exclusions = sum(workflow_excludes_zero), .groups = "drop")
write.csv(sensitivity, file.path(output, "module_sensitivity_summary.csv"), row.names = FALSE)
tex_table(data.frame(Setting = gsub("_", "\\\\_", sensitivity$specification), Lambda = sensitivity$lambda_gamma,
                    Bandwidth = sensitivity$h_multiplier, Workflows = sensitivity$workflows,
                    Intervals = sensitivity$intervals, Score = sensitivity$score_exclusions, Workflow = sensitivity$workflow_exclusions),
          c("Setting", "$\\lambda$", "$c_h$", "Fits", "Intervals", "Score", "Workflow"),
          "Completed module sensitivity workflows and interval exclusions of zero", "tab:joint-module-sensitivity", "module_sensitivity")
for (model in c("A", "B")) for (tau in c(.25, .5, .75)) {
  row <- reference[reference$model == model & reference$tau == tau & reference$module == "M4.11", ]
  prefix <- paste0("JointPlasma", model, c("Lower", "Median", "Upper")[match(tau, c(.25, .5, .75))])
  macro(paste0(prefix, "Estimate"), fmt(row$beta_hat))
  macro(paste0(prefix, "ScoreInterval"), interval(row$simultaneous_lower, row$simultaneous_upper))
  macro(paste0(prefix, "WorkflowInterval"), interval(row$workflow_simultaneous_lower, row$workflow_simultaneous_upper))
}
plot_data <- bind_rows(mutate(reference, interval_type = "Cluster score", lo = simultaneous_lower, hi = simultaneous_upper),
                       mutate(reference, interval_type = "Complete workflow", lo = workflow_simultaneous_lower, hi = workflow_simultaneous_upper))
plot_data$module <- factor(plot_data$module, levels = rev(c("M1.2", "M3.4", "M4.11", "M5.12", "M5.15")))
plot_data$model <- paste("Model", plot_data$model)
plot_data$quantile <- paste("Quantile", fmt(plot_data$tau, 2))
figure <- ggplot(plot_data, aes(x = module, y = beta_hat, ymin = lo, ymax = hi,
                               colour = interval_type, shape = interval_type)) +
  geom_hline(yintercept = 0, colour = "grey65", linewidth = .35) +
  geom_pointrange(position = position_dodge(width = .55), linewidth = .4, size = .3) +
  coord_flip() + facet_grid(model ~ quantile) +
  scale_colour_manual(values = c("Cluster score" = "#0072B2", "Complete workflow" = "#D55E00")) +
  scale_shape_manual(values = c("Cluster score" = 16, "Complete workflow" = 17)) +
  labs(x = NULL, y = "Profile coefficient for log C3", colour = NULL, shape = NULL) +
  theme_bw(base_size = 9, base_family = "sans") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        legend.position = "bottom", strip.background = element_rect(fill = "grey95"))
ggsave("manuscript/generated/figures/fig_joint_modules.pdf", figure, width = 7.1, height = 4.8, device = grDevices::cairo_pdf)
write.csv(plot_data, file.path(output, "module_figure_source.csv"), row.names = FALSE)

gene_files <- list.files("results/raw/application_exploratory/joint_profile_20261026", "^result[.]rds$", recursive = TRUE, full.names = TRUE)
if (length(gene_files) != 360L) stop("Incomplete independent gene roles.")
genes <- bind_rows(lapply(gene_files, function(path) readRDS(path)$table))
if (nrow(genes) != 1080L) stop("Incomplete role specific gene intervals.")
if (any(!is.finite(as.matrix(genes[c("beta_hat", "se", "lower", "upper")])))) stop("Unavailable exploratory interval values.")
macro("JointGeneRoles", length(gene_files))
macro("JointGeneIntervals", nrow(genes))
gene_summary <- genes %>% group_by(processing, model, tau, gene_id) %>% summarise(
  selected_roles = n(), coefficient_median = median(beta_hat),
  coefficient_q25 = quantile(beta_hat, .25), coefficient_q75 = quantile(beta_hat, .75), .groups = "drop") %>%
  arrange(processing, model, tau, desc(selected_roles), gene_id)
leading <- gene_summary %>% group_by(processing, model, tau) %>% slice_head(n = 3L) %>% ungroup()
write.csv(genes, file.path(output, "gene_role_intervals.csv"), row.names = FALSE)
write.csv(gene_summary, file.path(output, "gene_selection_summary.csv"), row.names = FALSE)
write.csv(leading, file.path(output, "gene_leading_summary.csv"), row.names = FALSE)
leading$route <- ifelse(leading$processing == "reference", "Reference", ifelse(grepl("official", leading$processing), "Raw arrays", "Median probes"))
for (model in c("A", "B")) {
  row <- leading[leading$processing == "reference" & leading$model == model & leading$tau == .25, ][1L, ]
  macro(paste0("JointLowerGene", model), row$gene_id)
  macro(paste0("JointLowerGeneCount", model), row$selected_roles)
}
tex_table(data.frame(Route = leading$route, Model = leading$model, Quantile = fmt(leading$tau, 2),
                    Gene = gsub("_", "\\\\_", leading$gene_id), Count = paste0(leading$selected_roles, "/20"),
                    Median = fmt(leading$coefficient_median), IQR = interval(leading$coefficient_q25, leading$coefficient_q75)),
          c("Processing", "Model", "$\\tau$", "Gene ID", "Selections", "Median", "Interquartile range"),
          "The three most frequently selected genes in each independent role experiment", "tab:joint-genes", "gene_selection", long = TRUE)
writeLines(unname(macros), file.path(generated, "inference_values.tex"))
jsonlite::write_json(list(method = "joint_profile_design_20261026", fixed_design = "joint_profile_fixed_design_20261026",
  generated = as.character(Sys.time()), independent_settings = nrow(registry), independent_datasets = nrow(full) / 4L,
  independent_intervals = nrow(full), paired_ablation_datasets = length(ablation_files),
  module_workflows = length(module_files), module_fits_including_deletions = length(module_files) * 130L,
  module_intervals = nrow(modules), exploratory_fits = length(gene_files), exploratory_intervals = nrow(genes),
  full_registered_simulation = "pending", source_tables = output,
  fixed_interface_checks = "results/scratch/joint_profile_manuscript_20261026/fixed_checks.csv",
  conic_solver_checks = "results/scratch/joint_profile_manuscript_20261026/conic_checks.csv",
  reproduction_checks = "results/scratch/score_validity_20261026/reproduction_checks.csv"),
  file.path(generated, "inference_manifest.json"), auto_unbox = TRUE, pretty = TRUE)
cat("Generated", nrow(paired), "paired simulation coordinates,", nrow(modules), "module intervals and", nrow(genes), "gene intervals\n")
print(overall)
print(reference[reference$score_excludes_zero | reference$workflow_excludes_zero,
                c("module", "model", "tau", "beta_hat", "simultaneous_lower", "simultaneous_upper", "workflow_simultaneous_lower", "workflow_simultaneous_upper")])
print(leading[leading$processing == "reference", c("model", "tau", "gene_id", "selected_roles")])
