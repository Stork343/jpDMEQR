import json
from pathlib import Path

import pandas as pd
import yaml

config = yaml.safe_load(Path("config/clinical_quantile_smooth_20261106.yml").read_text())
output = Path(config["report_directory"])
progress = pd.read_csv(output / "progress.csv")
verification = json.loads((output / "prediction_acceptance.json").read_text())
assert verification["status"] == "passed" and progress["status"].eq("completed").all()
comparison = pd.read_csv(output / "discovery_comparisons.csv")
numerical = pd.read_csv(output / "discovery_numerical_checks.csv")
reference = Path("results/preflight/clinical_quantile_20261105")
smooth_checks = json.loads((reference / "centered_smooth_acceptance.json").read_text())
centering = json.loads((reference / "centering_acceptance.json").read_text())
labels = {"ra_westchina_20261105": "RA", "pregnancy_detroit_20261105": "Pregnancy",
          "gnhs_baseline_bmi_20261105": "GNHS"}
decisions = []
for study_id, group in comparison[comparison["decision_quantile"]].groupby("study_id", sort=False):
    best = group.loc[group["comparator_loss"].idxmin()]
    old = group[group["comparator"].eq("JOINT_DIRECT")].iloc[0]
    decisions.append({"study_id": study_id, "tau": best["tau"], "best_comparator": best["comparator"],
                      "new_loss": best["current_loss"], "best_comparator_loss": best["comparator_loss"],
                      "relative_gain_over_best": best["relative_improvement"],
                      "relative_gain_over_joint_direct": old["relative_improvement"],
                      "linear_screen_passed": bool(best["relative_improvement"] >= 0.02),
                      "formal_validation_used": False, "manuscript_adoption": False})
decisions = pd.DataFrame(decisions)
decisions.to_csv(output / "screening_decisions.csv", index=False)
table = decisions[["study_id", "tau", "new_loss", "best_comparator", "best_comparator_loss", "relative_gain_over_best"]].copy()
table["study_id"] = table["study_id"].map(labels)
table["tau"] = table["tau"].map(lambda value: f"{value:.2f}")
for name in ["new_loss", "best_comparator_loss"]:
    table[name] = table[name].map(lambda value: f"{value:.6f}")
table["relative_gain_over_best"] = table["relative_gain_over_best"].map(lambda value: f"{value:.2%}")
table.columns = ["队列", "tau", "新方法损失", "既有最低损失方法", "该方法损失", "相对改善"]
pregnancy = decisions[decisions["study_id"].eq("pregnancy_detroit_20261105")].iloc[0]
timing = numerical[["study_id", "fits", "total_fitting_seconds"]].copy()
timing["study_id"] = timing["study_id"].map(labels)
timing["total_fitting_seconds"] = timing["total_fitting_seconds"].map(lambda value: f"{value / 60:.2f}")
timing.columns = ["队列", "完整拟合数量", "累计拟合分钟"]
inference_text = "患者正交推断正在计算。"
adaptive_path = output / "adaptive_inference.csv"
if adaptive_path.exists():
    adaptive = pd.read_csv(adaptive_path)
    display = adaptive[["study_id", "target", "estimate", "standard_error", "exchangeable_weight"]].copy()
    display["study_id"] = display["study_id"].map(labels)
    display.columns = ["队列", "检查系数", "估计值", "标准误", "交换型分量权重"]
    gnhs = adaptive[adaptive["study_id"].eq("gnhs_baseline_bmi_20261105")].iloc[0]
    ra = adaptive[adaptive["study_id"].eq("ra_westchina_20261105")].iloc[0]
    densities = pd.read_csv(output / "density_diagnostics.csv")
    ra_density = densities[densities["study_id"].eq("ra_westchina_20261105")]["ratio"]
    gnhs_density = densities[densities["study_id"].eq("gnhs_baseline_bmi_20261105")]["ratio"]
    inference_text = display.to_markdown(index=False, floatfmt=".6f") + f"""\n\nGNHS 的估计标准误由独立工作协方差下的 {gnhs['independence_standard_error']:.6f} 变为 {gnhs['standard_error']:.6f}，降低 {gnhs['estimated_se_reduction']:.2%}。这些系数用于方法计算检验，结果来自探索患者。

RA 中方差贡献最高的五位患者占估计方差的 {ra['largest_five_patient_variance_share']:.2%}，GNHS 对应比例为 {gnhs['largest_five_patient_variance_share']:.2%}。以训练分组外预测残差计算的零点密度与分位数间距密度作比较，后者平均值与前者之比在 RA 为 {ra_density.min():.2f} 至 {ra_density.max():.2f}，在 GNHS 为 {gnhs_density.min():.2f} 至 {gnhs_density.max():.2f}。该比值也随预测误差和核带宽变化。

独立重建全部20个投影方向的患者 influence 后，数值差异为零。"""
text = f"""# 临床分位数方法重设计与检验

## 估计目标

科学参数由临床结局的条件分位数模型定义。`CQ_SMOOTH` 使用中心化平滑损失估计初始回归，平滑程度为训练患者数的负五分之一次方。核的第一阶矩为零，因此在密度光滑和局部曲率条件下，平滑引起的总体近似误差为二阶。

全部预定临床变量进入可估计设计。蛋白质系数使用 elastic net，模型没有 floor(n/log n) 的入选变量上限。临床投影坐标和蛋白质系数同时加入随患者数增加而消失的数值 ridge，以确定唯一计算结果。

患者正交推断使用原始分位数指示得分。条件密度与工作协方差共同决定投影，标准误按完整患者计算。独立与交换型工作协方差的估计再按患者 influence 的协方差进行凸组合。

## 数学与数值检查

真实训练数据上的 {smooth_checks['checks']} 项独立二次规划检查全部通过，目标值最大差异为 {smooth_checks['maximum_objective_error']:.2e}，拟合预测的最大均方根差异为 {smooth_checks['maximum_prediction_rmse']:.2e}。{centering['cases']} 项确定性总体积分检查验证了核中心化及二阶近似误差，没有生成临床观测。

ORTHOGONAL_EXPANSION.md 给出患者交叉拟合的渐近展开，并列明初始预测与工具函数的误差乘积条件。THEORY_REVIEW.md 推导了正交工具函数的方差界。具体学习步骤的收敛速度和有限样本区间覆盖仍需进一步检验。

## 三个队列的完整比较

三个队列均已完成五个患者分组和三个分位数，共45项拟合。训练中的变量处理和参数选择按患者分组进行。每项拟合比较四个 elastic net 比例、四十个正则化参数，使用三个训练分组选择参数。

下表报告预定 decision_quantile，损失在每位患者的记录内平均，再对患者等权平均。表格依据已经完成的线性比较结果。三个预测分位数按既定规则排序并限制为非负。RA 的单位为 DAS28-CRP 评分，Pregnancy 为孕周，GNHS 为 BMI 的 kg/m²。

{table.to_markdown(index=False, disable_numparse=True)}

Pregnancy 相对旧 JOINT_DIRECT 的损失降低 {pregnancy['relative_gain_over_joint_direct']:.2%}。当前三个队列均未达到相对全部既有线性比较方法至少改善2%的探索条件。此次结果没有进入独立验证，manuscript_adoption=false。

## 计算时间与独立复算

{timing.to_markdown(index=False, disable_numparse=True)}

累计拟合时间为各个任务的耗时之和，三个进程并行计算。计时包含相应任务的参数选择与最终拟合，预先保存的患者分组和预处理参数按原登记使用。

R 根据保存系数独立重建全部预测，最大差异为 {verification['maximum_reconstruction_error']:.2e}。R 与 Python 的患者等权损失最大差异为 {verification['maximum_loss_difference']:.2e}。每个参数路径都检查原始坐标下的 KKT 条件。

## 患者正交推断

{inference_text}

## 文件与研究状态

完整模型位于 `{config['output_root']}`。本目录保存全部分位数损失、比较差异和数值核验结果。`screening_decisions.csv` 记录既定筛选条件的判定，`formal_validation_used` 与 `manuscript_adoption` 均为 false。

当前方法定义见 `docs/METHOD_SPECIFICATION.md` 和 `research/clinical_quantile_20261105/CENTERED_SMOOTHING.md`。
"""
(output / "REPORT_zh.md").write_text(text.strip() + "\n")
print(decisions.to_string(index=False))
