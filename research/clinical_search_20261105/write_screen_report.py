import argparse
import json
from pathlib import Path

import pandas as pd
import yaml

parser = argparse.ArgumentParser()
parser.add_argument("config", type=Path)
args = parser.parse_args()
config = yaml.safe_load(args.config.read_text())
output = Path(config["report_directory"])
status = json.loads((output / "discovery_status.json").read_text())
cohort = json.loads((output / "cohort_status.json").read_text())
counts = pd.read_csv(output / "cohort_counts.csv").set_index("role")
losses = pd.read_csv(output / "discovery_loss.csv")
wide = losses.query("scoring == 'primary'").pivot(index="method", columns="tau", values="loss")
ordered = wide.sort_values(config["decision_quantile"])
display = ordered.map(lambda x: f"{x:.4f}").reset_index()
display.columns = ["方法", "第25百分位", "第50百分位", "第75百分位"]
best = ordered.index[0]
target = status["best_screening_candidate"]
target_loss = wide.loc[target, config["decision_quantile"]]
best_loss = wide.loc[best, config["decision_quantile"]]
assert not status["continue_nonlinear_comparison"]
report = f"""# 类风湿关节炎蛋白质探索比较

## 数据与患者角色

West China Hospital 的公开临床表与蛋白质矩阵包含 {cohort['patients']} 名患者的 {cohort['records']} 次有效 DAS28-CRP 记录，{cohort['patients_with_repeats']} 名患者有重复测量，每名患者最多 {cohort['maximum_cluster_size']} 次。数据来源为 DOI {config['source_doi']}，采用已经包含2026年更正内容的公开文本。原始来源文件通过文件大小与压缩内容完整性检查。

患者角色在关联模型拟合前固定。探索部分有 {counts.loc['discovery', 'patients']} 名患者、{counts.loc['discovery', 'records']} 次记录，另外 {counts.loc['validation', 'patients']} 名患者、{counts.loc['validation', 'records']} 次记录保留验证。探索采用五个完整患者折，折内三个患者折选择参数。DAS28-CRP 的组成变量及 CRP 蛋白质均不进入输入。临床调整使用年龄、性别、ACPA 和随访阶段。

## 五折比较结果

全部 {status['completed_jobs']} 项计算完成。表中 pinball loss 对每名患者等权平均，数值越低表示预测误差越小。每个预测来自未包含该患者的训练样本；中位数填补和标准化参数均由相应训练患者确定。

{display.to_markdown(index=False, disable_numparse=True)}

预定决定分位数为0.75。{target} 的损失为 {target_loss:.4f}，{best} 为 {best_loss:.4f}，前者高出 {(target_loss / best_loss - 1) * 100:.2f}%。两种 JOINT 形式均未满足相对全部线性比较方法至少改善2%的探索条件。该候选按已登记条件结束探索，forest 与 XGBoost 比较没有启动。

## 研究状态

本次没有读取保留验证患者的模型结果，formal_validation_used 保持 false，正式验证次数没有增加。当前结果没有支持采纳这一应用，manuscript_adoption 保持 false。完整模型及交叉验证预测保存在 `results/raw/application/{config['study_id']}`，汇总值可由 `summarise_discovery.py` 重新生成。

原研究已经报告蛋白质与疾病活动度的线性关系及分组差异。本次没有开展蛋白质关联的正式检验，后续研究结论仍需要独立患者检验。
"""
(output / "REPORT_zh.md").write_text(report, encoding="utf-8")
print(output / "REPORT_zh.md")
