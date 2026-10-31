# 乳样实证分析结果

在纳入完整乳样光谱和牧场采样条件，并处理同一奶牛的重复测量后，这组数据未检出 A2 乳型与体细胞评分的明确关联。关于三个分位数之间的系数差异，整体 Wald 检验得到 p＝0.623，同样未达到显著标准。

分析使用 1,349 头奶牛的 2,253 份乳样，设计矩阵包含 1,073 个变量。研究对象是 A2 与 A1 乳型在这些可观测条件下的 profile 关联。18 项登记拟合均已完成，补充解释使用保存的奶牛影响值计算联合协方差。

## 调整后的乳型关联

三个分位数使用固定 `Lambda=1`。下列区间在这三个参数之间采用 Bonferroni 校正，置信水平为 95%。

| 分位数 | 估计值 | 标准误 | 同时置信区间 |
| --- | ---: | ---: | --- |
| 0.25 | −0.020 | 0.101 | [−0.263, 0.223] |
| 0.50 | 0.018 | 0.118 | [−0.264, 0.299] |
| 0.75 | 0.119 | 0.140 | [−0.216, 0.453] |

中位数系数为 0.018，点估计接近零，其区间仍允许 −0.264 至 0.299 个 SCS 单位的变化。现有估计的精度尚不足以确定这一关联的方向。

在条件模型之前，按奶牛等权计算的 A1 和 A2 乳样 SCS 中位数分别为 2.506 和 2.401。这些描述统计混合了不同的乳样光谱和采样条件，条件模型回答纳入这些信息后的乳型关联问题。

## 分位数之间的差异

三个分位数使用同一批奶牛，因此比较时保留了它们的估计相关性。全部成对差异如下，区间在三个比较之间采用 Bonferroni 校正。

| 系数差异 | 估计值 | 同时置信区间 |
| --- | ---: | --- |
| 0.50 减去 0.25 | 0.038 | [−0.241, 0.317] |
| 0.75 减去 0.25 | 0.139 | [−0.204, 0.482] |
| 0.75 减去 0.50 | 0.101 | [−0.241, 0.443] |

整体相等性检验的统计量为 0.947，自由度为 2，p＝0.623。上四分位数的点估计虽然较大，当前数据仍未检出其与下四分位数的明确差异。

## 敏感性结果

带宽乘数为 0.75 和 1.25 时，三个分位数系数相对参考设置的最大绝对变化约为 0.031 个 SCS 单位。工作惩罚为 0.5 和 2 的结果也均未排除零，各自对应其工作惩罚定义的 profile 参数。

限制为具有两次采样记录的 904 头奶牛后，中位数估计变为 −0.248，区间为 [−0.590, 0.094]。估计中心向负值移动，同时区间加宽。单次采样奶牛的平均 SCS 为 2.901，双次采样奶牛为 2.636，这项限制对应的奶牛总体具有不同的结局分布。全部 18 个区间均包含零，每种设置分别在其三个分位数之间进行校正。

全部估计和未取整的标准误保存在 `results/preflight/milk_scs_20261031/inference.csv`。

## 对统计方法的含义

这次应用完成了从完整光谱变量搜索到坐标区间的实际计算。相邻光谱变量的相关系数中位数为 0.9986，手稿的方法能够在本地计算机上处理这一规模的重复测量数据。推断输出给出了调整后乳型关联仍有多大不确定性，也显示了改变纳入奶牛范围时估计如何变化。

最终 18 项拟合队列约 3 分钟完成。完整坐标拟合时间的中位数为 17.0 秒，范围为 9.9 至 26.7 秒；共同使用的设计回归另外计算一次。各估计样本组包含 27 至 163 个选中变量，原始设计的全部 1,073 列参与变量搜索。

72 个带宽与样本组组合全部通过数值检查。原始坐标下最大 KKT 残差为 `5.17e-9`，最大 nuisance 梯度为 `9.99e-11`，QR 变换前后的最大线性预测差异为 `8.94e-9`。每个标准误均已从保存的奶牛影响值重新计算核验。

原始拟合文件及独立备份各有 18 份完整结果。工作结果目录为 `results/raw/application/milk_scs_20261031/covariance_qr_20261031_v2/`，备份目录为 `/Users/houjian/jpDMEQR-backups/milk_scs_20261031/covariance_qr_20261031_v2/`。

## 复现文件

分析方案见 `PROTOCOL_zh.md`，配置为 `config/application_milk_20261031.yml`。已保存结果可通过下列命令重新汇总并生成论文表格。

```sh
export TMPDIR="$PWD/results/scratch/runtime"
Rscript research/milk_application_20261031/summarise.R
Rscript research/milk_application_20261031/interpret_results.R
Rscript research/milk_application_20261031/build_tables.R
```

在具有原始数据且尚未创建本次结果目录的环境中，按以下顺序运行完整分析。

```sh
Rscript research/milk_application_20261031/prepare.R
Rscript research/milk_application_20261031/build_design_cache.R
Rscript research/milk_application_20261031/check_weighted_solver.R
Rscript research/milk_application_20261031/check_curvature_solver.R
Rscript research/milk_application_20261031/check_singleton_geometry.R
Rscript research/milk_application_20261031/run_remaining.R
Rscript research/milk_application_20261031/summarise.R
Rscript research/milk_application_20261031/interpret_results.R
Rscript research/milk_application_20261031/build_tables.R
```

手稿中的实证章节和 Supporting Information 均已使用这些结果。LaTeX 入口为 `manuscript/sjs_profile_quantile_main.tex` 和 `manuscript/sjs_profile_quantile_supplement.tex`，编译使用 `latexmk -xelatex`。
