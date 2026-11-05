import json
from pathlib import Path

import pandas as pd
import yaml

config = yaml.safe_load(Path("config/application_gnhs_baseline_bmi_20261105.yml").read_text())
output = Path(config["report_directory"])
status = json.loads((output / "discovery_status.json").read_text())
audit = json.loads((output / "screen_audit.json").read_text())
alignment = json.loads((output / "temporal_alignment_status.json").read_text())
cohort = json.loads((output / "cohort_status.json").read_text())
counts = pd.read_csv(output / "cohort_counts.csv").set_index("role")
assert audit["status"] == alignment["status"] == "passed"
assert status["completed_jobs"] == 75 and not status["formal_validation_used"]
risks = pd.read_csv(output / "discovery_loss.csv")
primary = risks.query("scoring == 'primary'").pivot(index="method", columns="tau", values="loss")
display = primary.sort_values(config["decision_quantile"]).map(lambda value: f"{value:.4f}").reset_index()
display.columns = ["方法", "第25百分位", "第50百分位", "第75百分位"]
comparison = pd.read_csv(output / "discovery_comparisons.csv")
decision = comparison.loc[comparison.scoring.eq("primary") & comparison.tau.eq(config["decision_quantile"])]
comparison_table = decision[["target", "comparator", "relative_gain"]].copy()
comparison_table.relative_gain = comparison_table.relative_gain.map(lambda value: f"{100 * value:.2f}%")
comparison_table.columns = ["研究方法", "比较方法", "相对改善"]
candidate = status["best_screening_candidate"]
gain = status["minimum_gain_by_candidate"][candidate]
if status["continue_nonlinear_comparison"]:
    decision_text = f"{candidate} 相对全部线性比较方法的最低改善为 {100 * gain:.2f}%，达到预先登记的2%探索条件。后续按照既定配置完成 quantile forest、XGBoost 和 LightGBM 比较。"
else:
    decision_text = f"两种 JOINT 形式均未达到预先登记的2%探索条件。{candidate} 的最低相对改善为 {100 * gain:.2f}%。该候选结束模型筛选，quantile forest、XGBoost 和 LightGBM 比较没有启动。"
report = f"""# 初次血清蛋白质与后续 BMI 预测

## 研究问题和数据

本研究使用 GNHS 初次检查的血清蛋白质和临床信息，预测 F2 与 F3 的 BMI。原始论文为 Cai 等的 Population serum proteomics uncovers a prognostic protein classifier for metabolic syndrome，DOI {config['source_doi']}。公开临床表中的 BMI 为连续测量，其他代谢指标以是否达到异常标准的标签提供。

匹配后包含 {cohort['patients']} 名受试者的 {cohort['records']} 条后续记录，每人最多两条。其中 {cohort['patients_with_repeats']} 人同时具有两个后续 BMI。探索队列有 {counts.loc['discovery', 'patients']} 人、{counts.loc['discovery', 'records']} 条记录，验证队列有 {counts.loc['validation', 'patients']} 人、{counts.loc['validation', 'records']} 条记录，公开患者编号没有重合。缺失结局没有填补。

两次后续 BMI 构成每人的结局向量。预测变量固定为初次检查的蛋白质、BMI、年龄、性别和代谢指标，并加入指定随访时长及预定时间乘积项。两个检测队列共同测得的413个蛋白质形成候选面板。这个设计采用结局向量之前确定的临床与蛋白质信息，对应手稿给定完整设计时的联合结局分布。

## 比较方法

探索队列采用按受试者划分的五折交叉验证，内层三个受试者分组选择参数。每个训练分组分别确定可观测比例筛选、填补与标准化参数。JOINT_DIRECT 直接拟合响应，JOINT_OFFSET 在 ridge 背景预测上拟合分位数修正。比较包含 exact quantile Lasso 及其 offset 形式、elastic net 和 ridge 的分位数校准，以及仅包含既有临床信息的 quantile regression。

预测统一限制为非负，并按三个分位数排序。pinball loss 在每人的记录内平均，再对所有受试者等权平均，单位为 kg/m²。第75百分位用于决定是否继续验证，第25和第50百分位共同报告。未排序预测作为敏感性结果保存。

## 探索结果

{display.to_markdown(index=False, disable_numparse=True)}

相对改善按比较方法损失减去研究方法损失，再除以比较方法损失计算。正数表示研究方法误差较低。

{comparison_table.to_markdown(index=False, disable_numparse=True)}

{decision_text}

## 数值复核与计算

所有75项拟合完成。从保存的模型独立重建预测，再由 R 重新计算受试者损失，与 Python 汇总的最大差异为 {audit['independent_loss_error']:.2e}。原始蛋白质及临床表的时间匹配也已独立核验，检查覆盖 {alignment['records_checked']} 条记录。

任务使用四个进程计算。preparation_summary.csv 保存各训练分组的预处理和 ridge 背景准备计时，discovery_runtime.csv 保存使用这些缓存参数后各个模型任务的计时。

## 研究状态

{counts.loc['validation', 'patients']} 名验证受试者的模型结果尚未读取，formal_validation_used 和 manuscript_adoption 均为 false。

数据来源及文件信息记录在 data_manifest.json。discovery_loss.csv 保存全部分位数损失，discovery_comparisons.csv 保存配对方法差异，temporal_alignment_checks.csv 和 prediction_reconstruction_checks.csv 保存数值核验结果。完整模型位于 `results/raw/application/gnhs_baseline_bmi_20261105`。
"""
(output / "REPORT_zh.md").write_text(report, encoding="utf-8")
print(output / "REPORT_zh.md")
