# 夜间血糖预测的临床数据检验

IOBP2 随机试验和 Loop 前瞻性研究均提供真实的纵向血糖记录。预测时点为 IOBP2 随机分组后第 14 日结束，或 Loop 入组后第 14 日结束。仅用此时已经取得的信息，预测随后三个 21 日时段的夜间高值。每夜取最高 10% 读数的均值，每个时段对合格夜间等权平均。第 75 百分位指这个连续结局在条件相近患者和时段中的较高水平。两项研究各自建立和评价模型。

## 预测比较

IOBP2 的探索患者中有 160 人符合监测要求，共提供 480 个结局。CQ_SMOOTH 的患者等权分位数预测损失为 7.4488。相比临床分位数回归的 7.6262，损失降低了 2.33%。经分位数校准的 elastic net 达到 7.1041，CQ_SMOOTH 的损失相对增加了 4.85%。此前登记的治疗前预测分别得到 28 名和 82 名监测完整的探索患者，未达到预定的 100 人要求。这里的 160 人回答的是治疗两周后的预测问题。

Loop 有 252 名合格探索患者，每人提供三个结局，共 756 个。CQ_SMOOTH 的损失为 6.7914，比校准 elastic net 的 7.3754 降低了 7.92%。临床分位数回归已经达到 6.8248，两者仅相差 0.0334，改善为 0.49%。研究登记要求相对每项比较方法均降低损失至少 2%。

## 完整曲线提供的信息

两周临床指标已经包含夜间高值、第二周的夜间高值及其相对第一周的变化。额外的 288 个变量描述一天中各时段的血糖分布。在两项研究各自的五组训练结果中，都有四组把这 288 个系数全部估为零。其余一组在 IOBP2 保留两个变量，在 Loop 保留一个。

补充比较在相同患者和相同训练组中，仅改变校准 elastic net 能否使用完整曲线。在 IOBP2，加入曲线变量使损失由 6.9969 升至 7.1041。在 Loop，损失由 7.0562 升至 7.3754。两项研究均未显示这些曲线变量能在已有临床指标的基础上改善该模型的预测。

## 稿件使用决定

两个公开数据包均已完成下载及压缩包校验，分析程序重建了全部探索结局和曲线变量。两项研究均满足数据结构与样本条件，各自的预测比较均未达到采纳标准。IOBP2 留出的 114 名患者及 Loop 留出的 216 名患者均只用于名册分组，其随访结局仍未用于正式验证。两项分析均作为方法开发记录保存，稿件尚无经独立验证的夜间血糖发现。逐项结果见 [IOBP2 计算报告](../clinical_quantile_iobp2_20261106/REPORT_zh.md) 和 [Loop 计算报告](../clinical_quantile_loop_20261106/REPORT_zh.md)。

IOBP2 的资料来自 Bionic Pancreas Research Group 发表的多中心随机试验，DOI 10.1056/NEJMoa2205225。Loop 的资料来自 Lum 等发表的前瞻性研究，DOI 10.1089/dia.2020.0535。开放数据由 Jaeb Center for Health Research 发布。

The source of the data is the Insulin Only Bionic Pancreas Pivotal Trial, but the analyses, content and conclusions presented herein are solely the responsibility of the authors and have not been reviewed or approved by the Bionic Pancreas Research Group or Beta Bionics.

The source of the data is the Loop Study (sponsored by the Jaeb Center for Health Research and funded by the Helmsley Charitable Trust), but the analyses, content and conclusions presented herein are solely the responsibility of the authors and have not been reviewed or approved by the study sponsor.
