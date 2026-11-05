# 固定参数网格上的路径计算

Quantile Lasso 在训练交叉验证中计算完整的 40 个参数。选中参数的位置确定之后，最终重拟合仍使用相同的完整参数网格，并按原来的顺序计算至选中位置。每个已经计算的位置保留全变量的 primal dual gap 和可行性检查。

maximum_solutions 只指定需要求解的网格位置数。完整网格的生成、初始模型与各位置的数值设置保持一致。最后一个已求解位置之后的约束清理不影响该位置已经保存的系数。

check_pooled_prefix.R 使用 IBDMDB 和 DrugMatrix 的真实训练记录，与已保存的完整路径比较系数、预测和最优性检查。IBDMDB 检查第 1、5、17、40 个位置，DrugMatrix 检查第 1 和第 5 个位置，并检查通用预测接口的选中参数映射。检查记录保存于 results/preflight/empirical_selection_20261101/pooled_prefix。
