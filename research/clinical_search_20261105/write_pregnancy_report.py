import json
from pathlib import Path

import pandas as pd
import yaml

config = yaml.safe_load(Path("config/application_pregnancy_detroit_20261105.yml").read_text())
output = Path(config["report_directory"])
status = json.loads((output / "discovery_status.json").read_text())
audit = json.loads((output / "screen_audit.json").read_text())
assert audit["status"] == "passed" and not status["continue_nonlinear_comparison"]
counts = pd.read_csv(output / "cohort_counts.csv").set_index("role")
loss = pd.read_csv(output / "discovery_loss.csv").query("scoring == 'primary'").pivot(index="method", columns="tau", values="loss")
table = loss.sort_values(0.5).map(lambda value: f"{value:.4f}").reset_index()
table.columns = ["方法", "第25百分位", "第50百分位", "第75百分位"]
target = status["best_screening_candidate"]
best = loss[0.5].idxmin()
report = f"""# 孕期蛋白质的孕周分位数探索

公开数据包含166名孕妇的666次采样，每名孕妇均有2至6次重复测量。数据来源为 DOI 10.1371/journal.pone.0181468，1125个蛋白质测量采用作者公开的供应商校准及极端值截取版本。结局为采样孕周，后续子痫前期诊断仅用于患者分层分组，模型输入只包含蛋白质。

探索部分包含 {counts.loc['discovery', 'patients']} 名孕妇的 {counts.loc['discovery', 'records']} 次采样，另有 {counts.loc['validation', 'patients']} 名孕妇的 {counts.loc['validation', 'records']} 次采样保留验证。五个完整患者折分别用于评价，折内三个患者折选择参数。判定分位数在拟合前固定为0.5。

## 预测比较

全部75项计算完成。下表为每名孕妇等权平均的 pinball loss，单位为周，较低数值表示误差较小。CLINICAL_QR 仅包含截距。

{table.to_markdown(index=False, disable_numparse=True)}

{target} 的第50百分位损失为 {loss.loc[target, 0.5]:.4f}，{best} 为 {loss.loc[best, 0.5]:.4f}，前者高出 {(loss.loc[target, 0.5] / loss.loc[best, 0.5] - 1) * 100:.2f}%。两种 JOINT 形式均未达到相对所有线性比较方法至少改善2%的探索条件。forest 与 XGBoost 比较没有启动，66名保留验证孕妇的模型结果尚未读取。

## 研究判断

已有研究已经使用 elastic net 构建孕期蛋白质时钟，也已经报告蛋白质丰度随孕周变化的10%、50%和90%分位数曲线。相应原文为 DOI 10.1016/j.ajog.2017.12.208 和 DOI 10.1016/j.ajog.2017.02.037。本次检验多个蛋白质条件下的孕周分位数估计，当前探索结果没有形成可供独立验证的方法优势。

从保存模型重建预测及独立计算患者损失的检查均通过，损失最大数值差异为 {audit['independent_loss_error']:.2e}。formal_validation_used 和 manuscript_adoption 均为 false，正式验证次数没有增加。完整模型及患者预测保存在 `results/raw/application/pregnancy_detroit_20261105`。
"""
(output / "REPORT_zh.md").write_text(report, encoding="utf-8")
print(output / "REPORT_zh.md")
