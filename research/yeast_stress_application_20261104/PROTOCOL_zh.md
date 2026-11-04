# 浓度系列的生长分布比较

本次登记使用 Peter 等人 2018 年公开数据中具有连续浓度或温度设置的培养条件。完整清单包含 anisomycin、caffeine、cycloheximide、formamide、benomyl 和 40°C 至 42°C 的温度条件。每个菌株在相应系列中提供两个或三个观测。全部系列沿用盐胁迫分析在关联拟合前保存的 848 个菌株及其角色，训练、探索评价和预留验证的菌株集合保持一致。

研究问题是环境强度增加时，生长分布较低部分与基因组特征的关联如何变化。各系列以 0.25 分位数的菌株平均 pinball loss 进行筛选；数据变换、基因组输入和低维调整沿用已经登记的方式。培养条件采用完整分类效应，各条件与 SNP、ORF 家族拷贝数和调整变量均具有交互项。这样可以表达三个浓度之间的非线性变化。所有方法都接收相同的设计变量。JOINT 的 nuisance 为菌株截距。

先登记完整计算及评价规则。JOINT 与经训练残差校准的均值模型完成所有系列的比较。JOINT 的探索评价损失更低时，计算同一系列的 exact quantile Lasso；同时优于这两类模型的系列进入 clustered quantile forest 及相关公开 RF 方法的比较。完整记录每项已尝试比较，未经计算的方法保持未完成状态。选择任何正式验证候选前，固定候选顺序、科学问题及比较方法。

原论文已开展 GWAS；2025 年预印本 Predictive models of the genetic bases underlying budding yeast fitness in multiple environments 使用同一资源建立 RF 等预测模型，并分析环境及基因型交互。新的结论需要具体的分位数差异和预留菌株中的方法优势。基因名称、一般的环境交互关系和单个较低的探索损失均不足以完成采纳。当前配置的 manuscript_adoption 为 false。
