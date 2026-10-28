# C3局部曲率审计规则

本规则于2026年10月22日在表达特征与C3关联拟合前确定，是本项目对既有局部曲率检查的计算定义。数值与判定方式由本轮项目统计审阅确定。

## 计算对象

对每个quantile与bandwidth组合，计算每名subject在该位置的平均Epanechnikov kernel weight，再对subject等权平均。quantile位置使用合格就诊的log C3经验quantile，经验quantile采用R的type 7。bandwidth为`h_multiplier * n^(-3/10)`，n为独立subject数量。

三个既定quantile为0.25、0.50、0.75，与0.75、1、1.25三个bandwidth multiplier组成九个组合。0.10敏感性quantile的三个bandwidth组合作为另一个集合，单独报告。

## Subject bootstrap

每次有放回抽取n个完整subject，保留每个subject的全部合格就诊。每次重抽样均重新计算log C3经验quantile和等权subject kernel均值。重复10000次，随机种子为20261022。

既定九个组合采用family confidence 0.95，分别取bootstrap均值分布的`0.05/9`分位数作为单侧percentile lower bound。敏感性三个组合使用`0.05/3`，经验分位数计算同样采用type 7。程序记录随机数算法、重复数量及每个组合的lower bound。

## 通过条件

九个既定组合的lower bound均需大于零。每个组合的local subject数量和Kish effective subject数量均需至少为`d_max+1`，其中`d_max=floor(n/log(max(n,3)))`是既有低维fixed-effect上限。在n为129时，该条件对应27。敏感性组合单独计算其通过状态。

Kish数量由subject平均kernel weight计算，即权重总和的平方除以权重平方和。该数量描述kernel信息在subject之间的分布。这个条件与边际C3分布有关，实际clinical/module设计继续接受已有reduced profile Hessian正曲率与数值条件检查。

## 统计前提

Bootstrap以独立subject为抽样单位，允许subject内部的就诊相关。其近似精度依赖subject抽样能够代表研究人群、局部kernel均值分布能够被subject重抽样近似，以及给定样本的quantile位置具有可重复性。由于quantile与kernel均值在每次bootstrap内共同计算，位置估计的变动进入bootstrap分布。Bonferroni分配覆盖的是本文件预先确定的组合。

该门禁提供边际C3局部信息的运行检查。profile设计的可识别性与区间校准由相应理论和数值检查处理。程序不使用表达特征与C3的关联结果选择这些数值。
