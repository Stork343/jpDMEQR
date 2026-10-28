#!/usr/bin/env Rscript

# 手稿表格与图形的构建脚本。输入为汇总目录 results/tables/<run_id>，
# 输出为 LaTeX 片段与 PDF 图形。任何数字都来自汇总表，不允许手工填入。
# 构建检查依据 docs/RESULTS_CONTRACT.md 第 9 节：每行绑定 experiment_id；
# 打印的重复次数等于冻结注册表的计数。

full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else
  "scripts/simulation/10_build_manuscript_assets.R"
root <- if (file.exists("R/metrics_v2.R")) "." else
  normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
source(file.path(root, "scripts", "00_source_v2.R"))
source_v2_module(root, "profile_v2", envir = environment())
source_v2_module(root, "simulation_v2", envir = environment())
source_v2_module(root, "metrics_v2", envir = environment())
source(file.path(root, "scripts", "simulation", "_run_registry_helpers.R"))

cli <- parse_cli_args_v2(commandArgs(trailingOnly = TRUE))
run_id <- cli$`run-id` %||% stop("Supply --run-id=<summary directory name>.")
registry_path <- cli$config %||% file.path(root, "config", "simulation_main_b200.csv")
mapping_path <- cli$mapping %||% file.path(root, "config", "manuscript_scenarios.csv")
out_dir <- cli$`out-dir` %||% file.path(root, "manuscript", "generated")
primary_method <- cli$`primary-method` %||% "PROFILE-DQR"
tables_dir <- file.path(root, "results", "tables", run_id)
if (!dir.exists(tables_dir)) stop("Missing summary directory: ", tables_dir)
figure_dir <- file.path(out_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

read_table <- function(name, required = TRUE) {
  path <- file.path(tables_dir, name)
  if (!file.exists(path)) {
    if (isTRUE(required)) stop("Missing summary table: ", path)
    return(data.frame())
  }
  read.csv(path, stringsAsFactors = FALSE)
}

write_fragment <- function(lines, name) {
  path <- file.path(out_dir, name)
  writeLines(lines, path, useBytes = TRUE)
  message("Wrote ", path)
}

tex_number <- function(x, digits = 3L) {
  out <- vapply(x, function(value) {
    if (!length(value) || !is.finite(value)) return("--")
    if (abs(value) < 1e-4) return(formatC(0, format = "f", digits = digits))
    formatC(value, format = "f", digits = digits)
  }, character(1))
  out
}

tex_sci <- function(x, digits = 2L) {
  vapply(x, function(value) {
    if (!length(value) || !is.finite(value)) return("--")
    formatC(value, format = "e", digits = digits)
  }, character(1))
}

tex_combined <- function(estimate, mcse, digits = 3L) {
  paste0(tex_number(estimate, digits), " (", tex_number(mcse, digits), ")")
}

tex_interval <- function(lower, upper, digits = 3L) {
  paste0("[", tex_number(lower, digits), ", ", tex_number(upper, digits), "]")
}

tex_label <- function(x) {
  gsub("_", "\\\\_", as.character(x))
}

# --- 输入与构建检查 ---------------------------------------------------------
registry <- read_simulation_registry_v2(registry_path)
mapping <- read.csv(mapping_path, stringsAsFactors = FALSE)
if (!all(c("row_label", "experiment_id") %in% names(mapping))) {
  stop("The scenario mapping needs row_label and experiment_id columns.")
}
unknown <- setdiff(mapping$experiment_id, registry$experiment_id)
if (length(unknown)) {
  stop("Scenario mapping refers to unknown experiment_id: ",
       paste(unknown, collapse = ", "))
}
accounting <- read_table("registry_accounting.csv")
if (!all(accounting$experiment_id %in% registry$experiment_id)) {
  stop("registry_accounting.csv contains rows outside the frozen registry.")
}
registry_total <- sum(as.integer(registry$replications) * vapply(
  registry$methods,
  function(x) length(parse_pipe_character_v2(x)),
  integer(1)
))
if (sum(accounting$scheduled_replications) != registry_total) {
  stop("Scheduled replications do not match the frozen registry count.")
}
if (any(accounting$missing_replications > 0L)) {
  stop("The summary directory does not cover every registered replication.")
}
cat(sprintf("Frozen registry replications: %d; recorded rows: %d\n",
            registry_total, sum(accounting$recorded_replications)))

method_summary <- read_table("method_summary.csv")
coverage_by_scenario <- read_table("coverage_by_scenario.csv")
coverage_table <- read_table("coverage_summary.csv")
inference_table <- read_table("inference_summary.csv")
uniformity <- read_table("null_pvalue_uniformity.csv")
studentised <- read_table("studentised_quantiles.csv")
selection_table <- read_table("selection_stability.csv")
theory_summary <- read_table("theory_summary.csv")
rate_slopes <- read_table("rate_slopes.csv", required = FALSE)
computation <- read_table("computation_summary.csv", required = FALSE)

scenario_meta <- merge(
  mapping,
  registry[, c("experiment_id", "module", "n_clusters", "p", "s", "tau", "q",
               "target_mode", "replications")],
  by = "experiment_id", all.x = TRUE, sort = FALSE
)
if (anyNA(scenario_meta$module)) stop("Scenario metadata is incomplete.")
scenario_meta <- scenario_meta[match(mapping$experiment_id, scenario_meta$experiment_id), ,
                               drop = FALSE]

lookup_row <- function(table, experiment_id, method_id) {
  hit <- table[table$experiment_id == experiment_id & table$method_id == method_id, ,
               drop = FALSE]
  if (nrow(hit) != 1L) {
    stop("Expected exactly one row for ", experiment_id, " / ", method_id,
         " but found ", nrow(hit), ".")
  }
  hit
}

# --- 表 1：推断结果汇总 -----------------------------------------------------
main_lines <- c(
  "\\begin{table}[t]",
  "\\caption{Primary inferential results for the proposed estimator. Coverage is reported at 95\\% and 90\\% with a replicate-level Monte Carlo standard error; averages are taken over the target coordinates of each scenario.}",
  "\\label{tab:simulation-results}",
  "\\small",
  "\\begin{tabularx}{\\textwidth}{@{}Xrrrrrrrr@{}}",
  "\\toprule",
  "Scenario & $\\tau$ & $n$ & Bias (MCSE) & Emp.\\ SD & Avg.\\ SE & SE/SD & Cov.\\ 95\\% (MCSE) & Cov.\\ 90\\% \\\\",
  "\\midrule"
)
for (ii in seq_len(nrow(scenario_meta))) {
  meta <- scenario_meta[ii, , drop = FALSE]
  cov_row <- lookup_row(coverage_by_scenario, meta$experiment_id, primary_method)
  main_lines <- c(main_lines, sprintf(
    "%s & %s & %d & %s & %s & %s & %s & %s & %s \\\\",
    tex_label(meta$row_label),
    tex_number(meta$tau, 2L),
    as.integer(meta$n_clusters),
    tex_combined(cov_row$bias, cov_row$bias_mcse),
    tex_number(cov_row$empirical_sd),
    tex_number(cov_row$mean_se),
    tex_number(cov_row$se_sd_ratio, 2L),
    tex_combined(cov_row$coverage_replication, cov_row$coverage_replication_mcse),
    tex_number(cov_row$coverage_replication_secondary)
  ))
}
main_lines <- c(main_lines,
  "\\bottomrule",
  "\\end{tabularx}",
  "\\end{table}"
)
write_fragment(main_lines, "tab_simulation_main.tex")

# --- 表 2：情景与方法对照 ---------------------------------------------------
method_lines <- c(
  "\\begin{table}[t]",
  "\\caption{Estimation and selection performance by scenario and method. Every row is bound to a frozen registry identifier; the final column reports successful replications out of the registered count.}",
  "\\label{tab:simulation-methods}",
  "\\small",
  "\\begin{tabularx}{\\textwidth}{@{}Xlrrrrrr@{}}",
  "\\toprule",
  "Scenario & Method & $\\ell_1$ (MCSE) & TPR & FDP & Pinball & Skill & Success \\\\",
  "\\midrule"
)
for (ii in seq_len(nrow(scenario_meta))) {
  meta <- scenario_meta[ii, , drop = FALSE]
  rows <- method_summary[method_summary$experiment_id == meta$experiment_id, ,
                         drop = FALSE]
  rows <- rows[order(match(rows$method_id, parse_pipe_character_v2(
    registry$methods[registry$experiment_id == meta$experiment_id][1]
  ))), , drop = FALSE]
  for (jj in seq_len(nrow(rows))) {
    method_lines <- c(method_lines, sprintf(
      "%s & %s & %s & %s & %s & %s & %s & %d/%d \\\\",
      if (jj == 1L) tex_label(meta$row_label) else "",
      tex_label(rows$method_id[jj]),
      tex_combined(rows$mean_l1[jj], rows$mcse_l1[jj]),
      tex_number(rows$mean_tpr[jj], 2L),
      tex_number(rows$mean_fdp[jj], 2L),
      tex_number(rows$mean_pinball[jj]),
      tex_number(rows$mean_pinball_skill_pooled[jj]),
      as.integer(rows$success[jj]),
      as.integer(rows$scheduled[jj])
    ))
  }
}
method_lines <- c(method_lines,
  "\\bottomrule",
  "\\end{tabularx}",
  "\\end{table}"
)
write_fragment(method_lines, "tab_simulation_methods.tex")

# --- 表 3：推断校准 ---------------------------------------------------------
inference_lines <- c(
  "\\begin{table}[t]",
  "\\caption{Inference calibration for the methods that construct an interval. Type I error is evaluated on null coordinates with a replicate-level Monte Carlo standard error; the bracket gives the binomial interval over pooled coordinates. Power is evaluated on active coordinates. Methods that produce point estimates and selections without an interval procedure are reported in the estimation table only.}",
  "\\label{tab:simulation-inference}",
  "\\small",
  "\\begin{tabularx}{\\textwidth}{@{}Xlrrrrr@{}}",
  "\\toprule",
  "Scenario & Method & Type I (MCSE) & Type I interval & Power & KS & Coverage \\\\",
  "\\midrule"
)
for (ii in seq_len(nrow(scenario_meta))) {
  meta <- scenario_meta[ii, , drop = FALSE]
  methods <- parse_pipe_character_v2(
    registry$methods[registry$experiment_id == meta$experiment_id][1]
  )
  # 只报告带区间程序的比较方法；没有区间程序的方法出现在估计与选择表中，
  # 缺失的原因写在构建日志中，不从表中静默删除。
  has_interval <- vapply(methods, function(m) {
    nrow(inference_table[inference_table$experiment_id == meta$experiment_id &
                           inference_table$method_id == m, , drop = FALSE]) == 1L
  }, logical(1))
  if (any(!has_interval)) {
    message("No interval procedure: ", meta$experiment_id, " / ",
            paste(methods[!has_interval], collapse = ", "))
  }
  methods <- methods[has_interval]
  if (!length(methods)) {
    stop("No method of ", meta$experiment_id, " carries an interval procedure.")
  }
  for (jj in seq_along(methods)) {
    inf_row <- lookup_row(inference_table, meta$experiment_id, methods[jj])
    ks_row <- uniformity[uniformity$experiment_id == meta$experiment_id &
                           uniformity$method_id == methods[jj], , drop = FALSE]
    cov_row <- coverage_by_scenario[coverage_by_scenario$experiment_id == meta$experiment_id &
                                      coverage_by_scenario$method_id == methods[jj], ,
                                    drop = FALSE]
    inference_lines <- c(inference_lines, sprintf(
      "%s & %s & %s & %s & %s & %s & %s \\\\",
      if (jj == 1L) tex_label(meta$row_label) else "",
      tex_label(methods[jj]),
      tex_combined(inf_row$type1_error_05, inf_row$type1_error_05_mcse),
      tex_interval(inf_row$type1_error_05_pooled_lower,
                   inf_row$type1_error_05_pooled_upper),
      tex_number(inf_row$power_05, 2L),
      tex_number(if (nrow(ks_row)) ks_row$ks_statistic[1] else NA_real_),
      tex_number(if (nrow(cov_row)) cov_row$coverage_replication[1] else NA_real_)
    ))
  }
}
inference_lines <- c(inference_lines,
  "\\bottomrule",
  "\\end{tabularx}",
  "\\end{table}"
)
write_fragment(inference_lines, "tab_simulation_inference.tex")

# --- 表 4：理论诊断 ---------------------------------------------------------
# 理论诊断依赖人口方向资产，因此覆盖注册表中带有 PROFILE-DQR-POP-H 的单元：
# 理论标度模块、profile 目标模块与精度压力模块。每一行都给出该单元的
# n、p、s 描述，数值全部来自 theory_summary.csv。
pop_h_rows <- registry[vapply(strsplit(registry$methods, "|", fixed = TRUE),
                              function(x) "PROFILE-DQR-POP-H" %in% x,
                              logical(1)), , drop = FALSE]
if (!nrow(pop_h_rows)) stop("No registry row carries PROFILE-DQR-POP-H.")
theory_lines <- c(
  "\\begin{table}[t]",
  "\\caption{Theorem diagnostics for the proposed estimator in the cells that carry a population effective direction. The target-approximation distance is reported for population-optimised targets; identity, nuisance-gradient and Dantzig entries are reported at working precision.}",
  "\\label{tab:simulation-theory}",
  "\\small",
  "\\begin{tabularx}{\\textwidth}{@{}Xrrrrrr@{}}",
  "\\toprule",
  "Cell ($n$, $p$, $s$) & Target distance & Schur identity & Nuisance gradient & Dantzig residual & Row error & Bahadur \\\\",
  "\\midrule"
)
theory_rows_written <- 0L
for (ii in seq_len(nrow(pop_h_rows))) {
  meta <- pop_h_rows[ii, , drop = FALSE]
  row <- theory_summary[theory_summary$experiment_id == meta$experiment_id &
                          theory_summary$method_id == primary_method, , drop = FALSE]
  if (nrow(row) != 1L) stop("Missing theory summary for ", meta$experiment_id, ".")
  theory_rows_written <- theory_rows_written + 1L
  theory_lines <- c(theory_lines, sprintf(
    "%s (%d, %d, %d) & %s & %s & %s & %s & %s & %s \\\\",
    tex_label(meta$experiment_id),
    as.integer(meta$n_clusters),
    as.integer(meta$p),
    as.integer(meta$s),
    tex_number(row$mean_beta_h_target_distance),
    tex_sci(row$max_profile_identity_error),
    tex_sci(row$max_nuisance_gradient),
    tex_sci(row$max_dantzig_residual),
    tex_sci(row$mean_precision_row_error),
    tex_number(row$mean_scaled_bahadur_remainder)
  ))
}
if (theory_rows_written != nrow(pop_h_rows)) {
  stop("The theory table did not cover every population-direction cell.")
}
theory_lines <- c(theory_lines,
  "\\bottomrule",
  "\\end{tabularx}",
  "\\end{table}"
)
write_fragment(theory_lines, "tab_simulation_theory.tex")

# --- 表 5：速率斜率 ---------------------------------------------------------
if (nrow(rate_slopes)) {
  slope_rows <- rate_slopes[rate_slopes$module == "A_theory_scaling", , drop = FALSE]
  if (nrow(slope_rows)) {
    slope_lines <- c(
      "\\begin{table}[t]",
      "\\caption{Empirical rate slopes from log-log regressions on the number of independent clusters, with the regression standard error.}",
      "\\label{tab:simulation-rates}",
      "\\small",
      "\\begin{tabularx}{\\textwidth}{@{}llrrr@{}}",
      "\\toprule",
      "Method & Metric & $p$ & Slope (MCSE) & Points \\\\",
      "\\midrule"
    )
    for (ii in seq_len(nrow(slope_rows))) {
      row <- slope_rows[ii, , drop = FALSE]
      slope_lines <- c(slope_lines, sprintf(
        "%s & %s & %d & %s & %d \\\\",
        tex_label(row$method_id),
        tex_label(row$metric),
        as.integer(row$p),
        tex_combined(row$slope, row$slope_mcse),
        as.integer(row$points)
      ))
    }
    slope_lines <- c(slope_lines, "\\bottomrule", "\\end{tabularx}", "\\end{table}")
    write_fragment(slope_lines, "tab_simulation_rates.tex")
  }
}

# --- 表 6：计算资源 ---------------------------------------------------------
if (nrow(computation)) {
  computation_lines <- c(
    "\\begin{table}[t]",
    "\\caption{Computation module. Timings are medians over the registered timed repeats; memory is the sampled resident set of the R process.}",
    "\\label{tab:simulation-computation}",
    "\\small",
    "\\begin{tabularx}{\\textwidth}{@{}Xlrrrrrr@{}}",
    "\\toprule",
    "Experiment & Method & Median time & IQR time & Max time & Peak memory (MB) & Iterations & Failures \\\\",
    "\\midrule"
  )
  for (ii in seq_len(nrow(computation))) {
    row <- computation[ii, , drop = FALSE]
    computation_lines <- c(computation_lines, sprintf(
      "%s & %s & %s & %s & %s & %s & %s & %s \\\\",
      tex_label(row$experiment_id),
      tex_label(row$method_id),
      tex_number(row$median_runtime, 1L),
      tex_number(row$q75_runtime - row$q25_runtime, 1L),
      tex_number(row$max_runtime, 1L),
      tex_number(row$max_peak_memory_mb, 1L),
      tex_number(row$median_outer_iterations, 0L),
      tex_number(row$failure_rate, 2L)
    ))
  }
  computation_lines <- c(computation_lines, "\\bottomrule", "\\end{tabularx}",
                         "\\end{table}")
  write_fragment(computation_lines, "tab_simulation_computation.tex")
}

# --- 表 7：注册表核算 -------------------------------------------------------
failure_lines <- c(
  "\\begin{table}[t]",
  "\\caption{Registry accounting. Every registered replication is accounted for; the flag marks scenarios with fewer than 95\\% successful replications.}",
  "\\label{tab:simulation-failures}",
  "\\small",
  "\\begin{tabularx}{\\textwidth}{@{}Xlrrrrrl@{}}",
  "\\toprule",
  "Experiment & Method & Registered & Recorded & Success & Failure & Not implemented & Flag \\\\",
  "\\midrule"
)
for (ii in seq_len(nrow(accounting))) {
  row <- accounting[ii, , drop = FALSE]
  failure_lines <- c(failure_lines, sprintf(
    "%s & %s & %d & %d & %d & %d & %d & %s \\\\",
    tex_label(row$experiment_id),
    tex_label(row$method_id),
    as.integer(row$scheduled_replications),
    as.integer(row$recorded_replications),
    as.integer(row$success),
    as.integer(row$failed),
    as.integer(row$not_implemented),
    if (isTRUE(row$low_success_flag)) "yes" else "no"
  ))
}
failure_lines <- c(failure_lines, "\\bottomrule", "\\end{tabularx}", "\\end{table}")
write_fragment(failure_lines, "tab_simulation_failures.tex")

# --- 表 8：选择稳定性 -------------------------------------------------------
selection_rows <- selection_table[selection_table$method_id == primary_method &
                                    selection_table$experiment_id %in% mapping$experiment_id, ,
                                  drop = FALSE]
if (nrow(selection_rows)) {
  selection_lines <- c(
    "\\begin{table}[t]",
    "\\caption{Selection stability across independent replications for the proposed estimator. Frequencies are the proportion of replications in which a coordinate is retained.}",
    "\\label{tab:simulation-selection}",
    "\\small",
    "\\begin{tabularx}{\\textwidth}{@{}Xlrrrrr@{}}",
    "\\toprule",
    "Scenario & Active freq. & Active $q_{90}$ & Stable active & Null freq. & Null $q_{90}$ & Unstable null \\\\",
    "\\midrule"
  )
  for (ii in seq_len(nrow(selection_rows))) {
    row <- selection_rows[ii, , drop = FALSE]
    label <- mapping$row_label[match(row$experiment_id, mapping$experiment_id)]
    selection_lines <- c(selection_lines, sprintf(
      "%s & %s & %s & %s & %s & %s & %s \\\\",
      tex_label(label),
      tex_number(row$active_selection_frequency_mean, 2L),
      tex_number(row$active_selection_frequency_q90, 2L),
      tex_number(row$stable_active_fraction, 2L),
      tex_number(row$null_selection_frequency_mean, 2L),
      tex_number(row$null_selection_frequency_q90, 2L),
      tex_number(row$unstable_null_count, 0L)
    ))
  }
  selection_lines <- c(selection_lines, "\\bottomrule", "\\end{tabularx}",
                       "\\end{table}")
  write_fragment(selection_lines, "tab_simulation_selection.tex")
}

# --- 图形 -------------------------------------------------------------------
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("ggplot2 is required for the manuscript figures.")
}
save_figure <- function(plot, name, width = 6.4, height = 4.0) {
  path <- file.path(figure_dir, name)
  ggplot2::ggsave(path, plot, width = width, height = height, device = "pdf")
  message("Wrote ", path)
}

module_a <- registry[registry$module == "A_theory_scaling", , drop = FALSE]
cov_a <- merge(
  module_a[, c("experiment_id", "n_clusters", "p", "s", "tau")],
  coverage_by_scenario[coverage_by_scenario$method_id == primary_method,
                       c("experiment_id", "coverage_replication",
                         "coverage_replication_mcse", "se_sd_ratio")],
  by = "experiment_id"
)
if (nrow(cov_a)) {
  cov_a$p_label <- factor(sprintf("p = %d", cov_a$p))
  coverage_plot <- ggplot2::ggplot(
    cov_a, ggplot2::aes(x = .data$n_clusters, y = .data$coverage_replication,
                        colour = .data$p_label)
  ) +
    ggplot2::geom_hline(yintercept = 0.95, linetype = 2, colour = "grey40") +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$coverage_replication - 1.96 * .data$coverage_replication_mcse,
                   ymax = .data$coverage_replication + 1.96 * .data$coverage_replication_mcse),
      width = 6
    ) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::labs(x = "Independent clusters", y = "Coverage at 95%",
                  colour = NULL) +
    ggplot2::theme_bw()
  save_figure(coverage_plot, "fig_coverage_clusters.pdf")

  ratio_plot <- ggplot2::ggplot(
    cov_a, ggplot2::aes(x = .data$n_clusters, y = .data$se_sd_ratio,
                        colour = .data$p_label)
  ) +
    ggplot2::geom_hline(yintercept = 1, linetype = 2, colour = "grey40") +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::labs(x = "Independent clusters",
                  y = "Average estimated SE / empirical SD", colour = NULL) +
    ggplot2::theme_bw()
  save_figure(ratio_plot, "fig_se_sd_ratio.pdf")
}

rate_a <- if (nrow(rate_slopes)) {
  rate_slopes[rate_slopes$module == "A_theory_scaling" &
                rate_slopes$metric %in% c("l1_error", "l2_error"), , drop = FALSE]
} else data.frame()
if (nrow(rate_a)) {
  rate_a$metric_label <- factor(rate_a$metric)
  rate_plot <- ggplot2::ggplot(
    rate_a, ggplot2::aes(x = .data$slope, y = .data$metric_label,
                         colour = .data$method_id)
  ) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = .data$slope - 1.96 * .data$slope_mcse,
                   xmax = .data$slope + 1.96 * .data$slope_mcse),
      height = 0.2
    ) +
    ggplot2::geom_point() +
    ggplot2::labs(x = "Empirical log-log slope in n", y = NULL, colour = NULL) +
    ggplot2::theme_bw()
  save_figure(rate_plot, "fig_rate_slopes.pdf")
}

if (nrow(inference_table)) {
  power_plot <- ggplot2::ggplot(
    inference_table,
    ggplot2::aes(x = .data$type1_error_05, y = .data$power_05,
                 colour = .data$method_id)
  ) +
    ggplot2::geom_hline(yintercept = 0.05, linetype = 2, colour = "grey40") +
    ggplot2::geom_vline(xintercept = 0.05, linetype = 2, colour = "grey40") +
    ggplot2::geom_point() +
    ggplot2::labs(x = "Type I error at 5%", y = "Power at 5%", colour = NULL) +
    ggplot2::theme_bw()
  save_figure(power_plot, "fig_type1_power.pdf")
}

uniform_rows <- uniformity[uniformity$experiment_id == mapping$experiment_id[1] &
                             uniformity$method_id == primary_method, , drop = FALSE]
if (nrow(uniform_rows)) {
  bin_columns <- grep("^bin_[0-9]+$", names(uniform_rows), value = TRUE)
  expected_columns <- grep("^expected_[0-9]+$", names(uniform_rows), value = TRUE)
  bin_data <- data.frame(
    lower = seq(0, 0.9, by = 0.1),
    fraction = as.numeric(uniform_rows[1, bin_columns]),
    expected = as.numeric(uniform_rows[1, expected_columns])
  )
  uniformity_plot <- ggplot2::ggplot(bin_data, ggplot2::aes(x = .data$lower)) +
    ggplot2::geom_col(ggplot2::aes(y = .data$fraction), width = 0.09,
                      fill = "grey70", colour = "grey30") +
    ggplot2::geom_hline(yintercept = 0.1, linetype = 2, colour = "black") +
    ggplot2::labs(x = "Null-coordinate p-value", y = "Fraction") +
    ggplot2::theme_bw()
  save_figure(uniformity_plot, "fig_null_pvalue_uniformity.pdf")
}

studentised_rows <- studentised[studentised$experiment_id == mapping$experiment_id[1] &
                                  studentised$method_id == primary_method, ,
                                drop = FALSE]
if (nrow(studentised_rows)) {
  quantile_plot <- ggplot2::ggplot(studentised_rows,
                                   ggplot2::aes(x = .data$normal, y = .data$empirical)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey40") +
    ggplot2::geom_point(size = 0.8) +
    ggplot2::labs(x = "Standard normal quantile",
                  y = "Studentised statistic quantile") +
    ggplot2::theme_bw()
  save_figure(quantile_plot, "fig_studentised_quantiles.pdf")
}

if (nrow(computation)) {
  time_plot <- ggplot2::ggplot(
    computation,
    ggplot2::aes(x = .data$median_outer_iterations, y = .data$median_runtime,
                 colour = .data$method_id)
  ) +
    ggplot2::geom_point() +
    ggplot2::labs(x = "Median outer iterations", y = "Median runtime (s)",
                  colour = NULL) +
    ggplot2::theme_bw()
  save_figure(time_plot, "fig_computation_time.pdf")
}

# 构建完整性：登记表要求写出的片段与图形必须全部存在，任何一项缺失就地中止，
# 避免在排版阶段才发现产物没有生成。
expected_assets <- expected_manuscript_assets_v2(registry)
missing_assets <- c(
  expected_assets$fragments[!file.exists(file.path(out_dir, expected_assets$fragments))],
  expected_assets$figures[!file.exists(file.path(figure_dir, expected_assets$figures))]
)
if (length(missing_assets)) {
  stop("Manuscript assets are incomplete: ", paste(missing_assets, collapse = ", "))
}
cat(sprintf("Manuscript asset completeness: %d fragments and %d figures present.\n",
            length(expected_assets$fragments), length(expected_assets$figures)))

message("Manuscript assets written to: ", out_dir)