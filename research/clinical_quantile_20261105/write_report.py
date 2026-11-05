import json
from pathlib import Path

import pandas as pd

root = Path("results/preflight/clinical_quantile_20261105")
progress = pd.read_csv(root / "progress.csv")
comparison = pd.read_csv(root / "discovery_comparisons.csv")
verification = json.loads((root / "prediction_acceptance.json").read_text())
engine = json.loads((root / "engine_acceptance.json").read_text())
inference = json.loads((root / "inference_implementation_acceptance.json").read_text())
assert verification["status"] == engine["status"] == inference["status"] == "passed"
labels = {"ra_westchina_20261105": "RA", "pregnancy_detroit_20261105": "Pregnancy",
          "gnhs_baseline_bmi_20261105": "GNHS"}
shown = comparison[comparison["decision_quantile"] & comparison["comparator"].isin([
    "JOINT_DIRECT", "ENET", "CLINICAL_QR"])].copy()
shown["study_id"] = shown["study_id"].map(labels)
shown["relative_improvement"] = shown["relative_improvement"].map(lambda value: f"{value:.2%}")
shown = shown[["study_id", "tau", "comparator", "current_loss", "comparator_loss", "relative_improvement"]]
shown.columns = ["队列", "tau", "比较方法", "CQ ENET 损失", "比较方法损失", "相对改善"]
display = progress[["study_id", "completed", "expected", "status"]].copy()
display["study_id"] = display["study_id"].map(labels)
display.columns = ["队列", "已完成拟合", "计划拟合", "状态"]
orthogonal = pd.read_csv(root / "ra_westchina_20261105_orthogonal_summary.csv")
orthogonal["method"] = orthogonal["method"].map({"independence": "独立工作协方差", "exchangeable": "交换型工作协方差"})
orthogonal.columns = ["计算方式", "患者数", "系数", "标准误", "区间下限", "区间上限"]
adaptive = pd.read_csv(root / "adaptive_inference.csv").set_index("study_id").loc["ra_westchina_20261105"]
text = f"""# 临床分位数方法开发记录

## 参数和估计方法

估计对象由临床结局的条件分位数模型定义。CQ ENET 最小化患者等权 pinball loss，对预定临床变量保留完整可估计设计，对蛋白质等其余变量使用 elastic net。参数选择采用既有患者分组，不设置 floor(n/log n) 的入选变量上限。

CQ ORTH 对预定系数计算患者正交得分。条件密度和工作协方差共同确定目标变量的加权投影，标准误由完整患者贡献计算。方法定义与条件写在 METHOD.md，推导检查写在 THEORY_REVIEW.md。条件分位数模型成立时，科学参数不随工作协方差改变。

## 实际计算

每个队列计划完成五个患者折和三个分位数，共十五项拟合。每项拟合使用三个训练分组选择参数，全部候选分组各包含四个 elastic net 比例和四十个正则化参数。

{display.to_markdown(index=False)}

点估计通过 {engine['objective_checks']} 项目标函数检查。正交得分的矩阵恒等式也已在 GNHS 的真实设计上核验；完整计算流程使用 1412 名训练受试者与 353 名评估受试者完成运行测试，该测试保留临床设计并使用十二个蛋白质检查接口。RA 随后完成全部候选变量下的五折正交计算。

完成队列的预测均由保存系数独立重建，最大差异为 {verification['maximum_reconstruction_error']:.2e}。R 与 Python 的患者等权损失最大差异为 {verification['maximum_loss_difference']:.2e}。

## 探索预测结果

下面使用原始登记的 decision_quantile。三个分位数统一排序并限制为非负，比较方法使用相同的评分规则。相对改善为比较方法损失减去 CQ ENET 损失，再除以比较方法损失。RA 的损失单位为 DAS28-CRP 评分，Pregnancy 为周。

{shown.to_markdown(index=False, floatfmt='.6f')}

Pregnancy 的结果衡量完整估计流程改变后的表现。CQ ENET 相对旧 JOINT_DIRECT 有所改善，ENET 的损失仍更低。RA 尚未取得预测改善。当前结果继续用于方法开发。

## 患者正交计算

RA 使用第75百分位和预定临床列 C_age 检查完整计算。两个工作协方差使用相同的初始分位数预测，按各自协方差重新计算投影。系数恢复为每一岁对应的 DAS28-CRP 评分差值，表中报告95%正态区间。

{orthogonal.to_markdown(index=False, floatfmt='.6f')}

交换型工作协方差在这次计算中的标准误更大。两个区间均包含零。这里报告探索队列的计算结果，区间的名义覆盖需要对应模型条件和专门的重复实验检验。

当前方法还根据两种患者 influence 的协方差计算凸组合，同时组合系数估计和标准误。在两个分量满足共同参数的中心化渐近展开、协方差估计一致时，组合的渐近方差不超过两个分量中较小的方差。RA 的交换型分量权重为 {adaptive['exchangeable_weight']:.3f}，组合标准误为 {adaptive['standard_error']:.6f}。该结果选择了独立工作协方差对应的估计。

## 理论与计算仍需完成的内容

总体分位数的识别和正交导数已经写出。高维 elastic net、条件密度间距和目标投影的具体收敛速度仍需联合推导，现有一般误差乘积条件尚不能替代这一部分。

GNHS 的完整参数选择耗时较长，计算正在进行。后续优化需要在相同目标和误差检查下完成。独立验证尚未启动，formal_validation_used=false，manuscript_adoption=false。

## 文件

方法规范为 docs/METHOD_SPECIFICATION.md，程序在 research/clinical_quantile_20261105。当前计算根目录为 results/raw/clinical_quantile_qp_20261105，聚合结果保存在本报告目录。原 profile 方法定义单独保存在 docs/methods/joint_profile_20261026.md。
"""
(root / "REPORT_zh.md").write_text(text.strip() + "\n")
print(f"报告已生成，完成 {progress['completed'].sum()} / {progress['expected'].sum()} 项拟合")
