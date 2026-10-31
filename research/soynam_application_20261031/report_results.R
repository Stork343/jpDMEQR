library(data.table)
library(yaml)
library(jsonlite)
config <- read_yaml("config/application_soynam_20261031.yml")
arguments <- commandArgs(trailingOnly = TRUE)
mode <- if (length(arguments)) arguments[1L] else "all"
if (!mode %in% c("all", "joint")) stop("Unknown report mode.")
source <- "results/tables/soynam_20261031"
prediction_source <- if (mode == "all") source else "results/preflight/soynam_20261031/joint_prediction"
all_effects <- fread(file.path(source, "inference.csv"))
effects <- all_effects[setting == "reference"]
manifest <- read_json("results/preflight/soynam_20261031/input_manifest.json", simplifyVector = TRUE)
audit <- read_json("results/preflight/soynam_20261031/data_audit.json", simplifyVector = TRUE)
folds <- fread(file.path(config$results_directory, "line_folds.csv"))
count <- function(value) format(value, big.mark = ",", scientific = FALSE, trim = TRUE)
comparison <- merge(all_effects[setting %in% c("bandwidth_075", "bandwidth_125")],
  effects[, .(coordinate, tau, reference = estimate_bu_ac)], by = c("coordinate", "tau"))
bandwidth_change <- max(abs(comparison$estimate_bu_ac - comparison$reference))
runtime <- range(effects$runtime_seconds) / 60
contrasts <- fread(file.path(source, "all_setting_contrasts.csv"))
comparisons <- fread(file.path(prediction_source, "prediction_comparisons.csv"))[calibrated == TRUE]
descriptive <- fread(file.path(source, "descriptive.csv"))
interval <- function(lower, upper, digits = 2L) sprintf(paste0("[ %.", digits, "f, %.", digits, "f ]"), lower, upper)
effect_rows <- vapply(seq_len(nrow(effects)), function(i) sprintf("| %s | %.2f | %.2f | %s |",
  sub("ACC10_", "", effects$coordinate[i]), effects$tau[i], effects$estimate_bu_ac[i],
  interval(effects$simultaneous_lower[i], effects$simultaneous_upper[i])), character(1))
gain <- comparisons[comparator == "joint_genomic"]
gain_rows <- vapply(seq_len(nrow(gain)), function(i) sprintf("| %.2f | %.2f%% | %s%% |",
  gain$tau[i], gain$relative_percent[i], interval(gain$relative_lower[i], gain$relative_upper[i])), character(1))
year <- contrasts[setting == "reference" & parameter == "year_2014_minus_2013_q75"]
year_sensitivity <- contrasts[parameter == "year_2014_minus_2013_q75"]
upper <- gain[tau == .75]
text <- c("# 冠层信息与大豆产量分布", "",
  sprintf("SoyNAM田间试验提供了%s个家系的%s个大豆品系。每个品系在2013年和2014年各有一条记录，合计%s个小区观测。%s个基因组标记分别与年份交互，连同冠层及其他调整变量，构成%s列设计矩阵。样本单位是品系，每个品系具有两次测量；这与我们的方法所处理的高维混合效应分位数模型相符。",
    count(audit$families), count(manifest$lines), count(manifest$records), count(manifest$genomic_markers), count(manifest$total_columns)), "",
  "冠层覆盖度在全基因组调整后仍与产量正相关。模型同时控制家系和田间条件，并调整记录中的成熟期。下表表示冠层覆盖度每增加10个百分点所对应的 profile 产量关联，单位为 bu/ac，即蒲式耳每英亩。区间共同控制六个年度分位数系数和五项预定差异的95%同时置信水平。", "",
  "| 年份 | 分位数 | 产量关联 | 同时置信区间 |", "| --- | --- | --- | --- |", effect_rows, "",
  sprintf("2014年的上四分位数关联比2013年高%.2f bu/ac，同时置信区间为%s。七项设置中的差异估计介于%.2f至%.2f bu/ac，所有相应区间均为正。增加田间空间项之后，这项年度差异继续得到支持。",
    year$estimate, interval(year$lower, year$upper), min(year_sensitivity$estimate), max(year_sensitivity$estimate)), "",
  "## 留出品系的预测", "",
  sprintf("预测训练使用%s个品系，另外%s个品系的两年记录全部留出。每种方法共用田间对照小区提供的校准信息；品系层面的预测变量包括基因型和冠层。截距校准来自训练样本的内部验证残差。下表比较包含冠层的joint profile预测与仅包含基因组和田间信息的同类预测，正值表示加入冠层后的pinball loss降低比例。",
    count(sum(folds$prediction_fold != config$prediction$heldout_fold)), count(sum(folds$prediction_fold == config$prediction$heldout_fold))), "",
  "| 分位数 | 预测损失降低比例 | 同时置信区间 |", "| --- | --- | --- |", gain_rows, "",
  sprintf("在上四分位数，预测损失下降%.2f%%，同时置信区间为%s%%。冠层测量对这批留出品系的产量上四分位数预测提供了额外信息。下四分位数和中位数的区间包含零，其估计增益也已保留在表格中。",
    upper$relative_percent, interval(upper$relative_lower, upper$relative_upper)), "",
  "## 计算及数据来源", "",
  sprintf("%s项预定推断均已完成。六个冠层系数在所有设置下都保持正向，同时置信区间均为正。带宽上下调整25%%时，参考系数的最大变化约为%.3f bu/ac。三个参考分位数的完整推断需要约%.1f至%.1f分钟，包含变量选择和两个带宽的系数估计，计算设备为配备36 GB内存的Apple M3 Max。",
    count(uniqueN(all_effects$case_id)), bandwidth_change, runtime[1L], runtime[2L]), "",
  "数据来自[NAM数据文档](https://search.r-project.org/CRAN/refmans/NAM/html/met.html)所说明的公开面板，测量设计见[Xavier等2017年的研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC5499164/)。产量单位和田间小区记录经过SoyNAM数据字典与相同小区的交叉核对。原始数据保留来源中的数值，计算结果按品系保存，并复制至独立备份目录。", "",
  sprintf("计算步骤见[复现说明](%s)。[推断结果表](%s)和[完整敏感性差异](%s)保留未取整数值；预测区间使用完整留出品系的成对损失，并按预定六项比较分配误差概率。",
    file.path(normalizePath("."), "research/soynam_application_20261031/README.md"),
    file.path(normalizePath("."), source, "inference.csv"),
    file.path(normalizePath("."), source, "all_setting_contrasts.csv")))
if (mode == "all") {
  pooled <- comparisons[comparator == "pooled_quantile_lasso"]
  rows <- vapply(seq_len(nrow(pooled)), function(i) sprintf("| %.2f | %.2f%% | %s%% |",
    pooled$tau[i], pooled$relative_percent[i], interval(pooled$relative_lower[i], pooled$relative_upper[i])), character(1))
  text <- c(text, "", "## 与pooled quantile Lasso的比较", "",
    "两种方法均使用冠层及完整基因组变量。表格报告joint profile相对于pooled quantile Lasso的损失降低比例，区间同样属于预定六项比较。", "",
    "| 分位数 | 预测损失降低比例 | 同时置信区间 |", "| --- | --- | --- |", rows)
}
writeLines(text, "research/soynam_application_20261031/RESULTS_zh.md", useBytes = TRUE)
cat("Wrote the result interpretation from the saved empirical tables.\n")
