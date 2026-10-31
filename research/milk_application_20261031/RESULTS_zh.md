# 乳样实证分析结果

2026 年 10 月 31 日完成全部 18 项登记分析。使用 1,349 头奶牛的 2,253 份乳样，完整设计矩阵包含 1,073 个变量。目标坐标是 A2 乳型与 A1 乳型的体细胞评分 profile 关联，调整变量包括全部 1,060 个 MIR 光谱变量及牧场采样条件。

## 参考设置

三个分位数使用固定 `Lambda=1`。下列区间在这三个参数之间采用 Bonferroni 校正，置信水平为 95%。

| 分位数 | 估计值 | 标准误 | 同时置信区间 |
| --- | ---: | ---: | --- |
| 0.25 | −0.020 | 0.101 | [−0.263, 0.223] |
| 0.50 | 0.018 | 0.118 | [−0.264, 0.299] |
| 0.75 | 0.119 | 0.140 | [−0.216, 0.453] |

三个区间都包含零，也包含一定幅度的正负关联。当前数据尚未确定调整后乳型关联的方向。这里的调整包含同期光谱所记录的乳样特征，系数按体细胞评分单位报告。

## 敏感性结果

五个附加设置覆盖了两个带宽、两个工作惩罚，以及具有两次采样记录的 904 头奶牛。每种设置分别对三个分位数进行 Bonferroni 校正，全部区间均包含零。

限制为双次采样奶牛后，中位数估计为 −0.248，区间为 [−0.590, 0.094]。这组数据对应具有两次采样记录的奶牛总体，其估计中心向负值移动。带宽设置的中位数估计分别为 0.049 和 0.016，区间宽度相近。不同工作惩罚分别定义各自的 profile 参数。

全部估计和未取整的标准误保存在 `results/preflight/milk_scs_20261031/inference.csv`。

## 数值运行

最终 18 项拟合队列约 3 分钟完成。完整坐标拟合时间的中位数为 17.0 秒，范围为 9.9 至 26.7 秒；共同使用的设计回归另外计算一次。各估计样本组包含 27 至 163 个选中变量，原始设计的全部 1,073 列参与变量搜索。

72 个带宽与样本组组合全部通过数值检查。原始坐标下最大 KKT 残差为 `5.17e-9`，最大 nuisance 梯度为 `9.99e-11`，QR 变换前后的最大线性预测差异为 `8.94e-9`。每个标准误均已从保存的奶牛影响值重新计算核验。

原始拟合文件及独立备份各有 18 份完整结果。工作结果目录为 `results/raw/application/milk_scs_20261031/covariance_qr_20261031_v2/`，备份目录为 `/Users/houjian/jpDMEQR-backups/milk_scs_20261031/covariance_qr_20261031_v2/`。

## 复现文件

分析方案见 `PROTOCOL_zh.md`，配置为 `config/application_milk_20261031.yml`。已保存结果可通过下列命令重新汇总并生成论文表格。

```sh
export TMPDIR="$PWD/results/scratch/runtime"
Rscript research/milk_application_20261031/summarise.R
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
Rscript research/milk_application_20261031/build_tables.R
```

手稿中的实证章节和 Supporting Information 均已使用这些结果。LaTeX 入口为 `manuscript/sjs_profile_quantile_main.tex` 和 `manuscript/sjs_profile_quantile_supplement.tex`，编译使用 `latexmk -xelatex`。
