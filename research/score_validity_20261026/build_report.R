source("research/score_validity_20261026/summarise_evaluation.R", local = TRUE)
if (run_id != "design_independent" || length(fits) != 1440L || nrow(table) != 5760L) stop("The independent evaluation is incomplete.")
if (length(failures)) stop("The independent evaluation contains failures.")
if (any(table$max_kkt > 1e-8) || any(!is.finite(as.matrix(table[, c("beta_tilde", "se", "lower", "upper")])))) stop("Invalid independent fit diagnostics.")
old_roots <- c("results/scratch/orthogonal_profile_score_20261025/independent_validation",
               "results/scratch/orthogonal_profile_score_20261025/dimensional_validation")
old_files <- unlist(lapply(old_roots, list.files, pattern = "^result[.]rds$", recursive = TRUE, full.names = TRUE))
old_table <- do.call(rbind, lapply(old_files, function(path) readRDS(path)$table))
old_table$type <- ifelse(old_table$truth == 0, "zero", "nonzero")
comparison <- do.call(rbind, lapply(c("nonzero", "zero"), function(type) {
  old <- old_table[old_table$type == type, ]
  new <- table[table$type == type, ]
  old_fraction <- aggregate(covered ~ case_id + replicate, old, mean)
  new_fraction <- aggregate(covered ~ case_id + replicate, new, mean)
  variance <- function(fractions) {
    groups <- split(fractions, fractions$case_id)
    total <- nrow(fractions)
    sum(vapply(groups, function(group) (nrow(group) / total)^2 * var(group$covered) / nrow(group), numeric(1)))
  }
  data.frame(type = type, old_coverage = mean(old$covered), new_coverage = mean(new$covered),
             difference = mean(new$covered) - mean(old$covered),
             difference_mcse = sqrt(variance(old_fraction) + variance(new_fraction)),
             old_mean_length = mean(old$upper - old$lower), new_mean_length = mean(new$upper - new$lower),
             length_ratio = mean(new$upper - new$lower) / mean(old$upper - old$lower))
}))
requests <- lapply(Sys.glob(file.path(root, "worker_request_*.rds")), readRDS)
completed <- lapply(Sys.glob(file.path(root, "worker_complete_*.rds")), readRDS)
if (length(requests) != 4L || length(completed) != 4L || sum(vapply(completed, `[[`, integer(1), "completed")) != 1440L) stop("Incomplete worker manifests.")
wall_minutes <- as.numeric(difftime(max(do.call(c, lapply(completed, `[[`, "finished"))),
                                     min(do.call(c, lapply(requests, `[[`, "started"))), units = "mins"))
target <- read.csv(file.path(directory, "selected_target_audit.csv"), stringsAsFactors = FALSE)
geometry <- read.csv("results/scratch/score_validity_20261026/weighted_profile_checks.csv")
interface <- read.csv("results/scratch/score_validity_20261026/design_joint_checks.csv")
if (any(geometry$error > geometry$tolerance) || any(interface$error > interface$tolerance)) stop("Numerical checks failed.")
validation <- data.frame(check_group = c("weighted_profile", "joint_interface"),
                          passed = c(nrow(geometry), nrow(interface)))
write.csv(comparison, file.path(directory, "comparison.csv"), row.names = FALSE)
write.csv(validation, file.path(directory, "validation_counts.csv"), row.names = FALSE)
write.csv(data.frame(datasets = 1440L, intervals = 5760L, workers = 4L,
                     wall_minutes = wall_minutes, max_kkt = max(table$max_kkt),
                     maximum_absolute_selected_target_shift = max(abs(target$target_shift)),
                     observed_target_shift_count = sum(target$material_shift)),
           file.path(directory, "run_metrics.csv"), row.names = FALSE)
scaling_path <- "results/scratch/score_validity_20261026/design_scaling/summary.csv"
if (!file.exists(scaling_path)) stop("The scaling experiment is incomplete.")
scaling <- read.csv(scaling_path)
if (nrow(scaling) != 4L || any(scaling$replications != 20L)) stop("Incomplete scaling replications.")
write.csv(scaling, file.path(directory, "scaling.csv"), row.names = FALSE)
format_percent <- function(x) sprintf("%.2f%%", 100 * x)
line <- function(data) paste0("| ", ifelse(data$type == "nonzero", "非零系数", "零系数"), " | ",
                              format_percent(data$coverage), " | ", format_percent(data$mcse), " | ",
                              format_percent(data$mc_lower), " 至 ", format_percent(data$mc_upper), " |")
scale_line <- function(data) sprintf("| %s | %.3f | %.3f | %.3f |", format(data$p, big.mark = ","), data$median, data$p90, data$maximum)
lines <- c(
  "# 独立验证结果", "",
  "本报告由 `build_report.R` 根据原始结果生成。方法标识为 `joint_profile_design_20261026`，入口函数为 `ds_fit_profile`。", "",
  "## 覆盖率", "",
  "24种设置各运行60次，共1,440个独立数据集和5,760个区间。设置包含200和400个独立样本组，变量数为500、1,000或2,000，检查三个分位数和两种随机效应结构。误差分布包含正态和标准化t3，另有较小回归系数及较高协变量相关性的设置。名义置信水平为95%。", "",
  "| 系数类型 | 覆盖率 | 蒙特卡洛标准误 | 覆盖率的蒙特卡洛95%区间 |",
  "| --- | ---: | ---: | ---: |",
  vapply(seq_len(nrow(overall)), function(i) line(overall[i, ]), character(1)), "",
  "同一数据集的两个坐标可能相关，汇总标准误按照数据集计算，并保留不同设置的分层结构。逐坐标的二项区间保存在 `coordinate_summary.csv`。", "",
  sprintf("1,440次计算全部完成，没有失败。四个R进程的实际运行时间为%.2f分钟，最终拟合的最大梯度残差为%.3g。", wall_minutes, max(table$max_kkt)), "",
  "## 计算时间", "",
  "以下结果在同一台本地计算机上使用一个R进程测量。每种变量数重复20次，均为200个独立样本组、随机截距、0.5分位数和四个目标系数。计时包含变量选择、估计和区间计算，代码编译与数据生成在计时范围之外。", "",
  "| 变量数 | 时间中位数 秒 | 时间90%分位数 秒 | 最长时间 秒 |",
  "| --- | ---: | ---: | ---: |",
  vapply(seq_len(nrow(scaling)), function(i) scale_line(scaling[i, ]), character(1)), "",
  "## 目标参数与区间长度", "",
  sprintf("根据各设置声明的高斯协变量分布，逐一计算了所选模型的总体目标。11,520个模型坐标中，目标变化超过1e-8的数量为%d，最大绝对变化为%.3g。生成系数只用于这一结果检查。", sum(target$material_shift), max(abs(target$target_shift))), "",
  sprintf("相同24种设置下，既有独立结果的非零系数覆盖率为%s，本次为%s。对应的平均区间长度从%.4f变为%.4f，增加%s。覆盖率高于名义水平的同时，区间长度也是需要报告的性能代价。", format_percent(comparison$old_coverage[1L]), format_percent(comparison$new_coverage[1L]), comparison$old_mean_length[1L], comparison$new_mean_length[1L], format_percent(comparison$length_ratio[1L] - 1)), "",
  "这两次比较采用相同设置和不同随机种子。逐设置结果及区间长度保存在 `scenario_summary.csv` 和 `coordinate_summary.csv`。", "",
  "## 数值检查", "",
  sprintf("%d项剖面导数和权重一致性检查通过；%d项联合估计接口检查通过，包含带宽协方差、变量单位变化和观测顺序变化。", nrow(geometry), nrow(interface)), "",
  "## 复现文件", "",
  "完整参数表为 `evaluation.csv`。独立数据种子由 `1800000000 + 10000 * case_index + replicate` 确定。原始结果保存在 `results/scratch/score_validity_20261026/design_independent`，汇总表保存在 `results/tables/score_validity_20261026/design_independent`。方法和推导分别见 `METHOD.md` 与 `JOINT_THEORY.md`。"
)
writeLines(lines, "research/score_validity_20261026/RESULTS.md", useBytes = TRUE)
cat("Independent wall minutes", wall_minutes, "\n")
print(comparison)
print(scaling)
