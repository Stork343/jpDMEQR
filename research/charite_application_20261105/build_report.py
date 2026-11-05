import json
import platform
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

config = yaml.safe_load(Path("config/application_charite_hb_20261105.yml").read_text())
output = Path(config["report_directory"])
raw = Path(config["output_root"])
source = Path(config["source_directory"])
download = json.loads((source / "manifest.json").read_text())
assert download["status"] == "downloaded_with_si" and len(download["files"]) == 13
audits = {name: json.loads((output / f"{name}_audit_status.json").read_text()) for name in ["python", "r"]}
assert all(item["status"] == "passed" for item in audits.values())
summary = pd.read_csv(output / "prediction_summary.csv")
comparisons = pd.read_csv(output / "prediction_comparisons.csv")
contrasts = pd.read_csv(output / "profile_contrasts.csv")
coefficients = pd.read_csv(output / "profile_coefficients.csv")
calibration = pd.read_csv(output / "quantile_calibration.csv")
jackknife = pd.read_csv(output / "jackknife_contrasts.csv")
counts = pd.read_csv(output / "cohort_counts.csv").set_index("role")
cohort = json.loads((output / "cohort_status.json").read_text())
inference = json.loads((output / "inference_status.json").read_text())
diagnostics = pd.read_csv(output / "inference_diagnostics.csv")
runtime = pd.read_csv(output / "training_runtime.csv")
crossing = pd.read_csv(output / "quantile_crossings.csv").set_index("method")
targets = pd.read_csv(output / "frozen_inference_targets.csv")
primary = summary.query("scoring == 'primary'").copy()
decision = comparisons.query("scoring == 'primary' and tau == 0.25").copy()
reference = contrasts.query("multiplier == 1").copy()
prediction_pass = bool(decision.superiority_pass.all())
association_pass = bool(reference.contrast_pass.any())
adoption = prediction_pass and association_pass
assert not adoption


def table(frame):
    return frame.to_markdown(index=False, disable_numparse=True)


def interval(lower, upper, percent=False):
    factor = 100 if percent else 1
    digits = 1 if percent else 3
    return f"[{factor * lower:.{digits}f}, {factor * upper:.{digits}f}]"


method_names = {
    "JOINT_OFFSET": "JOINT_OFFSET 研究方法",
    "JOINT_DIRECT": "JOINT_DIRECT",
    "POOLED": "quantile Lasso",
    "POOLED_OFFSET": "quantile Lasso offset",
    "ENET": "elastic net 分位数校准",
    "PRESS_RIDGE": "ridge 分位数校准",
    "CLINICAL_QR": "clinical quantile regression",
    "QRF_RAW": "quantile forest",
    "QRF_OFFSET": "quantile forest offset",
    "XGB_RAW": "XGBoost quantile regression",
    "XGB_OFFSET": "XGBoost quantile regression offset",
}
loss_matrix = primary.pivot(index="method", columns="tau", values="loss")
order = loss_matrix.sort_values(0.25).index
loss_table = pd.DataFrame({"方法": [method_names[x] for x in order]})
for tau in config["quantiles"]:
    loss_table[f"第{int(tau * 100)}百分位"] = [f"{loss_matrix.loc[x, tau]:.4f}" for x in order]
comparison_table = pd.DataFrame({
    "比较方法": decision.comparator.map(method_names),
    "相对改善百分比": decision.relative_gain.map(lambda x: f"{100 * x:.1f}"),
    "描述性95%区间": [interval(a, b, True) for a, b in zip(decision.lower95, decision.upper95)],
    "单侧98.75%下界": decision.one_sided_lower.map(lambda x: f"{100 * x:.1f}"),
})
contrast_table = pd.DataFrame({
    "蛋白质": reference.protein,
    "分位数对比": reference.contrast.map({"tau025_minus_tau050": "0.25 − 0.50", "tau025_minus_tau075": "0.25 − 0.75"}),
    "估计差异": reference.estimate.map(lambda x: f"{x:.3f}"),
    "多重比较区间": [interval(a, b) for a, b in zip(reference.lower, reference.upper)],
    "Bonferroni p值": reference.p_bonferroni.map(lambda x: f"{x:.3f}"),
})
runtime_table = runtime.groupby("method", sort=False).agg(拟合任务=("method", "size"), 秒数=("runtime_sec", "sum")).reset_index()
runtime_table["方法"] = runtime_table.method.map(method_names)
runtime_table["秒数"] = runtime_table["秒数"].map(lambda x: f"{x:.3f}")
runtime_table = runtime_table[["方法", "拟合任务", "秒数"]]
best = loss_matrix[0.25].idxmin()
target_loss = loss_matrix.loc[config["target_method"], 0.25]
best_loss = loss_matrix.loc[best, 0.25]
relative_excess = target_loss / best_loss - 1
joint_crossing = crossing.loc[config["target_method"]]
raw_decision = comparisons.query("scoring == 'raw' and tau == 0.25")
max_gradient = diagnostics[["full_kkt", "half_kkt", "full_nuisance_gradient", "half_nuisance_gradient"]].to_numpy().max()
python_error = max(audits["python"]["maximum_errors"].values())
inference_checks = pd.read_csv(output / "inference_reconstruction_checks.csv")
inference_error = inference_checks[["influence_error", "coefficient_error", "contrast_error", "se_error"]].to_numpy().max()
jackknife_pass = int(((jackknife.lower_jackknife > 0) | (jackknife.upper_jackknife < 0)).sum())
curve = calibration.query("scoring == 'primary' and method == 'JOINT_OFFSET'").sort_values("tau")
source_text = source / "PMC8201874.2.txt"
assert source_text.exists()

report = f"""# Charité 纵向血浆蛋白与血红蛋白分析

## 数据和研究问题

数据来自 Demichev 等发表于 Cell Systems 的 A time-resolved proteomic and prognostic map of COVID-19，DOI 10.1016/j.cels.2021.05.005。PMC version 2 提供的正文及全部补充文件共 {len(download['files'])} 个，已按官方文件大小核验。XLSX 和 ZIP 文件通过压缩文件完整性检查，两个 PDF 均可解析。完整来源清单见同目录 data_manifest.json。

Table S8 提供蛋白质与临床测量值，Table S9 给出患者和采样信息，Table S10 记录肽段数量。原始矩阵的 {cohort['source_rows']} 个采样编号均能匹配到患者资料。血红蛋白具有 {cohort['patients']} 名患者的 {cohort['records']} 次有效记录，单个患者最多 {cohort['maximum_cluster_size']} 次。结局保留 g/dL 单位，全部有效正值记录纳入分析，缺失结局没有填补。

研究检验两个问题。蛋白质与临床变量结合后，JOINT_OFFSET 能否改善新患者血红蛋白第25百分位的预测。其次，在固定临床条件与蛋白质设计下，蛋白质关联是否随血红蛋白分位数变化。第50和第75百分位同时报告。

原研究在采样当日缺少临床读数时使用此前两天以内的最近读数，公开附表没有每项检验的独立时间戳。因此，当前结局是与蛋白质采样配对的近期血红蛋白读数。数据没有提供本分析所需的逐次输血记录及患者性别字段，这两项因素无法纳入调整。

## 患者分组和比较方法

患者分组在查看模型关联结果之前完成。按首次有效记录的 WHO grade 是否达到6分层，固定随机种子20261105，训练部分包含 {counts.loc['training', 'patients']} 名患者、{counts.loc['training', 'records']} 次记录；验证部分包含 {counts.loc['validation', 'patients']} 名患者、{counts.loc['validation', 'records']} 次记录。训练与验证之间没有重复患者。三个完整患者折用于选择参数。

蛋白质候选要求至少三个肽段。移除含血红蛋白基因的 {cohort['excluded_response_protein_groups']} 个蛋白质组后，共有 {cohort['protein_candidates']} 个候选。蛋白质转换为 log2，训练患者加权可观测比例至少为80%的变量进入相应模型。中位数填补、缺失指示变量及标准化参数均由对应训练患者确定。临床变量包含年龄、BMI、Charlson index、入院后时间和当前 WHO grade。模型没有使用随后最高病情、患者结局或未来恶化状态。

JOINT_OFFSET 在 ridge 背景预测上估计分位数修正，JOINT_DIRECT 直接拟合响应。比较方法采用 exact quantile Lasso、elastic net、ridge、clinical quantile regression、grf quantile forest 和原生 XGBoost quantile regression，按协议加入相应 offset 版本。相同患者分组用于全部模型。原研究预测机械通气和 WHO grade，本分析将 XGBoost 用于连续血红蛋白的 quantile loss，当前结果评价这一共同预测任务。

预测值统一限制为非负，并将三个分位数排序。pinball loss 在每个患者内取平均，再对患者等权平均。数值越低表示预测误差越小。配对患者 bootstrap 重抽样20000次。预定判据要求第25百分位对全部10项比较方法的单侧98.75%改善区间下界均大于零。

## 预留患者的预测结果

{table(loss_table)}

JOINT_OFFSET 的第25百分位损失为 {target_loss:.4f}，{method_names[best]} 为 {best_loss:.4f}，前者高出 {100 * relative_excess:.1f}%。JOINT_OFFSET 对全部10项比较方法的点估计均没有改善，预定预测判据未通过。第50百分位的最低损失来自 {method_names[loss_matrix[0.5].idxmin()]}，第75百分位最低损失来自 {method_names[loss_matrix[0.75].idxmin()]}。

相对改善按“比较方法损失减去 JOINT_OFFSET 损失，再除以比较方法损失”计算，下表负数表示 JOINT_OFFSET 的损失更高。95%区间用于描述不确定性，正式判据使用最后一列。

{table(comparison_table)}

原始预测保留为敏感性分析，第25百分位仍有 {int((raw_decision.relative_gain < 0).sum())} 项比较的点估计没有改善。JOINT_OFFSET 原始预测在 {int(joint_crossing.crossing_records)}/{int(joint_crossing.records)} 次记录中出现分位数交叉，比例为 {100 * joint_crossing.crossing_records / joint_crossing.records:.1f}%。统一排序后，观察值低于第25、第50及第75百分位预测的患者加权比例分别为 {curve.fraction_below.iloc[0]:.3f}、{curve.fraction_below.iloc[1]:.3f} 和 {curve.fraction_below.iloc[2]:.3f}。较低分位数的覆盖比例接近目标值，pinball loss 仍高于比较方法，反映误差幅度也需要纳入评价。

## 蛋白质分位数关联

HP、HPX 和 TF 在协议中预先指定。训练数据根据第25与第75百分位的完整预测斜率差异选择 GSN、FBLN1 和 C6。全部目标在验证开始前固定。验证模型同时调整这六个蛋白质和临床变量，共 {inference['dimension']} 个可识别参数。

三个分位数共用 Lambda 等于1，响应与蛋白质尺度来自训练数据。估计采用 jp_fit_fixed 的两个带宽配对系数，标准误由患者 influence 及跨分位数协方差计算。默认带宽为 n^(-3/10)。推断对象为固定设计下的 regularised profile association，训练预测模型的完整斜率包含更广的蛋白质变量，二者的条件变量集合不同。

下表单位为蛋白质相对丰度每增加一倍时，血红蛋白关联差异的 g/dL 数值。{len(reference)}项对比共同使用错误概率{config['validation']['inference_alpha']}，区间按 Bonferroni 调整，逐项置信水平为{100 * (1 - config['validation']['inference_alpha'] / len(reference)):.4f}%。

{table(contrast_table)}

默认带宽的12项区间全部包含零。HPX 的第25与第75百分位差异具有最小的调整后 p 值，为 {reference.p_bonferroni.min():.3f}，仍高于预定阈值0.0125。40名验证患者提供的区间较宽，当前分析无法确立具有重复检验控制的分位数关联差异。

0.75倍和1.25倍带宽下共有 {int(contrasts.query('multiplier != 1').contrast_pass.sum())} 项对比通过对应标准。固定变换参数和带宽后逐个删除患者，jackknife 区间有 {jackknife_pass} 项排除零。敏感性结果与默认分析的判定一致。

## 与已有研究的关系

原论文已经报告血红蛋白、红细胞和红细胞压积与多项急性期蛋白相关，并讨论 HP、HPX 在 RRT 和 ECMO 患者中的降低。GSN 与 COVID 19 严重程度的关系以及 FBLN1 与年龄的关系也已出现在原研究。本次针对血红蛋白条件分位数的差异检验回答了不同问题，但全部正式检验均未通过，现有结果没有确立新增临床发现。

血红蛋白本身可通过常规检验获得。当前研究用途是比较估计方法和检验蛋白质关联；本分析未评价检测成本或临床治疗决策。

## 数值复核和计算时间

全部19项训练任务完成，保存了训练交叉验证结果和所选模型。R 记录的推断梯度及 KKT 残差最大为 {max_gradient:.2e}。独立 Python 计算重建了 {audits['python']['prediction_rows_checked']} 行原始及排序预测记录的损失、患者均值和相对改善，最大数值差异为 {python_error:.2e}。从保存的完整模型及患者得分重新计算 influence 和分位数对比，最大差异为 {inference_error:.2e}。

下表汇总各方法训练进程记录的秒数，包含参数选择和完整训练拟合。多个任务同时执行，表内秒数相加表示任务耗时之和，实际等待时间由任务启动安排决定。数据下载、预处理和验证推断时间另计。

{table(runtime_table)}

本次正式验证已完成，formal_validation_used 为 true。预测和蛋白质差异判据均未通过，manuscript_adoption 保持 false。研究代码、汇总表和本报告保存到 GitHub；患者记录、原始附表和完整模型保存在项目目录及独立备份目录。

## 结果文件

prediction_summary.csv 保存全部预测损失，prediction_comparisons.csv 保存配对区间。profile_coefficients.csv 包含各分位数系数及描述性95%区间，profile_contrasts.csv 保存正式对比。jackknife_contrasts.csv 和 inference_diagnostics.csv 提供敏感性结果及求解器记录。audit_status.json 汇集本次复核状态，study_status.json 记录研究判定。

数据来源全文保存在 data/raw/covid_charite_20261105/official_cloud/PMC8201874.2.txt，附表保存在同一目录。原文相关位置为 Results、Markers of RRT and ECMO、Charité patient cohort and clinical data、Prediction of current mechanical ventilation 及 Supplemental Information。
"""
(output / "REPORT_zh.md").write_text(report, encoding="utf-8")
(output / "data_manifest.json").write_text(json.dumps(download, indent=2) + "\n")
(output / "audit_status.json").write_text(json.dumps(audits, indent=2) + "\n")
status = {"study_id": config["study_id"], "status": "completed", "data_complete": True,
          "completed_training_jobs": len(runtime), "methods": primary.method.nunique(),
          "training_patients": int(counts.loc["training", "patients"]),
          "validation_patients": int(counts.loc["validation", "patients"]),
          "prediction_superiority_pass": prediction_pass, "association_contrast_pass": association_pass,
          "formal_validation_used": True, "manuscript_adoption": adoption,
          "completed_utc": datetime.now(timezone.utc).isoformat()}
(output / "study_status.json").write_text(json.dumps(status, indent=2) + "\n")
environment = {"python": platform.python_version(), "packages": {
    package: version(package) for package in ["numpy", "pandas", "PyYAML", "openpyxl", "pypdf", "tabulate", "statsmodels", "xgboost"]}}
(output / "python_environment.json").write_text(json.dumps(environment, indent=2) + "\n")
print(json.dumps(status, indent=2))
