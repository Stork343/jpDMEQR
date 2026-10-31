# 药物剂量反应数据研究

当前候选为 GDSC 药物筛选数据。研究问题是，基因表达与不同剂量下的反应不足有什么关系，分位数方法能否改善对较高存活读数的预测，以及发现的关联能否在另一轮筛选实验中重复。

GDSC1 具有多个剂量的连续测量，基因表达来自同一细胞系的 RNA-seq。实验单位依据 Cell Model Passports 的 `patient_id` 组成完整样本组。来自同一位患者的不同细胞系及衍生细胞系使用相同样本角色。角色划分使用元数据，全部关联计算开始之前保存，探索部分占患者编号的百分之七十。

GDSC 的孔位读数由同一实验板的对照孔进行标准化。采用 `gdscIC50::normalizeData` 的 `trim = FALSE`，保留实际标准化测量值。重复记录依照实验板及孔位标识核查。质量检查失败的孔位按来源程序处理。基因表达转换与变量选择分别使用对应训练样本。

每个来源与药物组合按实验日期选取最早的完整剂量曲线，具有五个或九个剂量点。这样每个来源对应同一实验周期。其余实验周期用于检验测量重复性，全部记录保存在数据对象中。探索训练、探索评价和正式验证遵循同一患者分组。

Koukouli 等人的 functional regression 使用基因与剂量的 B-spline 乘积，并通过变量筛选和 group SCAD 估计剂量变化的平均关联。其论文已经研究 BRAF 药物以及 MAPK 相关通路，原作者程序为 `fbioSelectr`。本研究需要与该程序进行共同样本上的比较，并纳入直接估计分位数的高维方法。评价使用完整患者组，每个方法使用相同原始信息，保留各自的训练过程。平均响应方法的分位数预测转换只使用训练样本。

统计方法的参数仍是固定 `Lambda` 的 profile 系数。剂量处的关联由相应系数及其联合协方差计算。分位数之间的差异进行直接检验，发现样本与验证样本分别保存。细胞系实验中的关联需要对应具体药物与浓度；研究结论不扩展至患者用药效果。

候选需要通过实际预测比较，并形成经文献核查和保留样本验证的科学结论，之后才可采纳。GDSC1 与 GDSC2 的测量技术分别记录；使用 iDRR 时读取研究标识完整的原始剂量记录。跨研究重复测量保留共同来源关系，验证角色在全部研究中一致。

参考来源为 [Koukouli 等人的研究](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008066)、[GDSC 实验方法](https://depmap.sanger.ac.uk/documentation/datasets/drug-sensitivity/)、[细胞系来源关系](https://depmap.sanger.ac.uk/documentation/cell-models/model-relationships/)和 [RNA-seq 处理说明](https://depmap.sanger.ac.uk/documentation/datasets/expression/)。
