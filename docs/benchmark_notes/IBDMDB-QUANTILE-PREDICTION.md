# 患者代谢组的分位数预测

数据来自 Lloyd-Price 等人的 IBDMDB 研究，原文使用混合模型分析代谢物与疾病活动，并利用调整后的残差研究多组学关联。公开资料包含每次采样的患者编号和 calprotectin 读数。本文按同一次采样匹配，预测新的完整患者的 log calprotectin。来源为 [Multi-omics of the gut microbial ecosystem in inflammatory bowel diseases](https://pmc.ncbi.nlm.nih.gov/articles/PMC6650278/)，DOI 10.1038/s41586-019-1237-9。

当前文件匹配到 85 名参与者的 349 次采样。17 名参与者保留用于正式验证，其余 68 名参与者的 269 次采样用于五折探索评价。按 pooled QC 变异系数选择非重复测量后，519 个有名称的代谢物符合预定质量条件。患者内的不同采样保持在同一个评价折，变量变换和背景学习均使用相应训练部分。

calprotectin 的两个公开字段在 97 次采样中可以核对，最大数值差异小于 0.0001。fecalcal_ng_ml 提供了单位名称和更完整的小数精度。代码在该字段可用时采用它，否则保留 fecalcal 的原始读数。预测损失计算使用自然对数尺度。

五折比较包含均值背景加 JOINT 修正、相同背景和变量下的 exact quantile Lasso、直接 JOINT 模型，以及经过分位数校准的 ridge 与 elastic net。第 75 分位数为登记的评价分位数。第 50 分位数的探索结果随后接受单独登记的完整比较，其余结果共同保留。

非线性比较直接调用 grf 2.6.1 的 quantile_forest。每个森林使用 2000 棵树，参数覆盖两种 splitting、两种 mtry、两种叶节点大小，以及 honesty 的两个设置。每个训练折从无校准、训练残差校准和 OOB 残差校准中选择。原始响应与均值残差各拟合一组森林。实现使用原生 clusters 和 equalize.cluster.weights 参数，文档为 [quantile_forest](https://grf-labs.github.io/grf/reference/quantile_forest.html) 和 [predict.quantile_forest](https://grf-labs.github.io/grf/reference/predict.quantile_forest.html)。

全部 55 项拟合完成后，报告脚本重新计算保存模型的预测、调参损失及选中校准量的加权分位数条件。原始模型与复算结果一致。探索性区间通过重新抽取完整患者计算，交叉验证训练部分之间有重叠，因此这些区间用于描述当前探索结果。

第 50 分位数的 JOINT 修正预测损失为 0.3865046，原始响应 quantile forest 为 0.3754642。当前候选未达到全部比较方法的优势条件，正式验证状态和手稿采纳状态均为 false。完整数值来自 results/preflight/empirical_selection_20261101/ibdmdb_forests。
