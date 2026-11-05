# 初次血清蛋白质与后续 BMI 预测

## 研究问题和数据

本研究使用 GNHS 初次检查的血清蛋白质和临床信息，预测 F2 与 F3 的 BMI。原始论文为 Cai 等的 Population serum proteomics uncovers a prognostic protein classifier for metabolic syndrome，DOI 10.1016/j.xcrm.2023.101172。公开临床表中的 BMI 为连续测量，其他代谢指标以是否达到异常标准的标签提供。

匹配后包含 2649 名受试者的 4524 条后续记录，每人最多两条。其中 1875 人同时具有两个后续 BMI。探索队列有 1765 人、3253 条记录，验证队列有 884 人、1271 条记录，公开患者编号没有重合。缺失结局没有填补。

两次后续 BMI 构成每人的结局向量。预测变量固定为初次检查的蛋白质、BMI、年龄、性别和代谢指标，并加入指定随访时长及预定时间乘积项。两个检测队列共同测得的413个蛋白质形成候选面板。这个设计采用结局向量之前确定的临床与蛋白质信息，对应手稿给定完整设计时的联合结局分布。

## 比较方法

探索队列采用按受试者划分的五折交叉验证，内层三个受试者分组选择参数。每个训练分组分别确定可观测比例筛选、填补与标准化参数。JOINT_DIRECT 直接拟合响应，JOINT_OFFSET 在 ridge 背景预测上拟合分位数修正。比较包含 exact quantile Lasso 及其 offset 形式、elastic net 和 ridge 的分位数校准，以及仅包含既有临床信息的 quantile regression。

预测统一限制为非负，并按三个分位数排序。pinball loss 在每人的记录内平均，再对所有受试者等权平均，单位为 kg/m²。第75百分位用于决定是否继续验证，第25和第50百分位共同报告。未排序预测作为敏感性结果保存。

## 探索结果

| 方法          | 第25百分位   | 第50百分位   | 第75百分位   |
|:--------------|:-------------|:-------------|:-------------|
| CLINICAL_QR   | 0.4352       | 0.5514       | 0.4537       |
| JOINT_OFFSET  | 0.4376       | 0.5530       | 0.4553       |
| JOINT_DIRECT  | 0.4378       | 0.5526       | 0.4567       |
| POOLED        | 0.4380       | 0.5589       | 0.4629       |
| ENET          | 0.4367       | 0.5578       | 0.4645       |
| POOLED_OFFSET | 0.4535       | 0.5644       | 0.4702       |
| PRESS_RIDGE   | 0.5367       | 0.6872       | 0.5741       |

相对改善按比较方法损失减去研究方法损失，再除以比较方法损失计算。正数表示研究方法误差较低。

| 研究方法     | 比较方法      | 相对改善   |
|:-------------|:--------------|:-----------|
| JOINT_DIRECT | POOLED        | 1.35%      |
| JOINT_DIRECT | POOLED_OFFSET | 2.87%      |
| JOINT_DIRECT | ENET          | 1.68%      |
| JOINT_DIRECT | PRESS_RIDGE   | 20.46%     |
| JOINT_DIRECT | CLINICAL_QR   | -0.66%     |
| JOINT_OFFSET | POOLED        | 1.66%      |
| JOINT_OFFSET | POOLED_OFFSET | 3.17%      |
| JOINT_OFFSET | ENET          | 1.99%      |
| JOINT_OFFSET | PRESS_RIDGE   | 20.71%     |
| JOINT_OFFSET | CLINICAL_QR   | -0.35%     |

两种 JOINT 形式均未达到预先登记的2%探索条件。JOINT_OFFSET 的最低相对改善为 -0.35%。该候选结束模型筛选，quantile forest、XGBoost 和 LightGBM 比较没有启动。

## 数值复核与计算

所有75项拟合完成。从保存的模型独立重建预测，再由 R 重新计算受试者损失，与 Python 汇总的最大差异为 1.28e-15。原始蛋白质及临床表的时间匹配也已独立核验，检查覆盖 4524 条记录。

任务使用四个进程计算。preparation_summary.csv 保存各训练分组的预处理和 ridge 背景准备计时，discovery_runtime.csv 保存使用这些缓存参数后各个模型任务的计时。

## 研究状态

884 名验证受试者的模型结果尚未读取，formal_validation_used 和 manuscript_adoption 均为 false。

数据来源及文件信息记录在 data_manifest.json。discovery_loss.csv 保存全部分位数损失，discovery_comparisons.csv 保存配对方法差异，temporal_alignment_checks.csv 和 prediction_reconstruction_checks.csv 保存数值核验结果。完整模型位于 `results/raw/application/gnhs_baseline_bmi_20261105`。
