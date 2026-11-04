# 酵母环境系列的相关研究

Peter 等人公开了 1011 个菌株的基因组信息及其中 971 个菌株的环境生长比例。论文使用 FaST-LMM 对变换后的表型开展关联分析。补充表 15 报告 ENA5、VAR1589 和 MDM36 与 1 M NaCl 条件的关联，1.5 M 条件没有通过其阈值的关联记录。原文及补充表已保存于 data/raw/yeast_1011_20261104。来源为 [Genome evolution across 1011 Saccharomyces cerevisiae isolates](https://www.nature.com/articles/s41586-018-0030-5)，DOI 10.1038/s41586-018-0030-5。

Žun 和 Petrovič 使用连续回交与等位基因替换实验研究盐耐受性，并利用原公开菌株资源分析 ASG1 与 ENA 剂量。论文在 1 M NaCl 条件下比较了不同 ASG1 类型和 ENA 剂量的菌株。新的统计分析以跨浓度、跨分位数的 profile 关联为研究对象。来源为 [Discovery and validation of ASG1 as a novel determinant of NaCl tolerance](https://pmc.ncbi.nlm.nih.gov/articles/PMC12774596/)，DOI 10.1093/g3journal/jkaf254。

Segura Abá 等人在 2026 年分析了同一资源中的 750 个二倍体菌株。其预测方法包括 rrBLUP、Bayesian LASSO、BayesC、RF 和 XGBoost。RF 使用训练部分的变量重要程度进行筛选；解释部分研究了不同环境和菌株的变量贡献，并报告 benomyl 相关基因交互。因此，新的候选需比较直接分位数方法及非线性预测方法，并估计具有具体尺度的分位数关联差异。来源为 [Predictive models of the genetic bases underlying budding yeast fitness in multiple environments](https://pmc.ncbi.nlm.nih.gov/articles/PMC13458404/)，DOI 10.1093/nargab/lqag084。

该研究的 RF 参数网格包含 max_depth 为 3、5、10，max_features 为 0.1、0.5、sqrt、log2、None，以及 n_estimators 为 100、500、1000。作者的 ML-Pipeline 使用 scikit-learn。项目比较采用完整菌株分组，变量筛选与残差校准均需要在对应训练部分完成。作者测试集、全队列生成的特征排名和已拟合模型的分数不用于本项目的预留验证。

当前模型采用共同的 SNP 和 ORF 家族拷贝数输入，按训练菌株计算变换与 PC。探索评价与坐标推断分别保存。Ridge 和 elastic net 提供经过培养条件校准的均值预测；exact quantile Lasso 与 clustered quantile forest 直接拟合分位数。RF 的方法来源、参数范围和变量选择程序需要在其比较结果中完整记录。
