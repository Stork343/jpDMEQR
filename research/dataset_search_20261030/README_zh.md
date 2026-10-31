# 奶牛乳样数据与高维分位数推断

## 数据及研究问题

本项目使用 Universitat Autònoma de Barcelona 公开的奶牛乳样数据，研究在给定完整红外光谱和牧场采样条件后，A2 乳型与体细胞评分分布的关联。目标是在不同分位数估计同一个预先指定的乳型系数，并计算其区间。数据同时保存奶牛编号、重复采样记录及连续数值结局。全部实证结果见 `research/milk_application_20261031/RESULTS_zh.md`。

数据集名称为 *Genetic protein fractions, chemical composition, and mid-infrared spectra of cow milk*，DOI 为 `10.34810/data2338`。原始 Excel 已从官方仓库完整下载，文件为 21,859,858 bytes，读取后得到 2,270 行观测和 1,073 列。公开说明记录了测量仪器：1,060 个 MIR 光谱变量来自 MilkoScan 7RM，体细胞计数来自 Fossomatic。[官方数据说明](https://portalrecerca.uab.cat/en/datasets/genetic-protein-fractions-chemical-composition-and-mid-infrared-s/)

体细胞评分与乳蛋白基因型之间已有可参照的研究问题。Čítek 等使用重复测量混合模型研究这一关联，并考察牧场及泌乳阶段等因素。这里拟研究给定高维光谱后的分位数 profile 系数，其定义包括同期乳样特征的调整。[Čítek 等 2022](https://doi.org/10.5713/ab.20.0720)

## 实际记录核查

核查使用原始 Excel 的数值和编号，结果由 `audit_milk.R` 生成。

| 记录范围 | 奶牛数量 | 乳样数量 | 每头采样次数 |
| --- | ---: | ---: | --- |
| 官方文件 | 1,356 | 2,270 | 1 至 3 |
| 编号及基因型一致且数据完整 | 1,349 | 2,253 | 1 至 2 |
| 具有两次完整采样 | 904 | 1,808 | 2 |

有 7 个奶牛编号涉及基因型变化、跨牧场记录或同次采样出现多条记录，合计 17 份乳样。这些编号的全部记录留在原始文件，建模清单排除相应奶牛。处理规则仅使用编号和元数据。445 头只有一次采样的奶牛保留在全体可用记录中，904 头双次采样奶牛另有清单。

光谱列为 `240` 至 `1299`，每列均有变化，全部数值有限。体细胞计数及所需元数据均无缺失，数值记录的单位是 cells/µL。转换采用 `SCS = log2(SCC/100) + 3`，并保持原始记录精度。可用数据有 572 个不同的 SCS 数值，出现最多的数值占 3.51%。Excel 中的采样日期覆盖 2023 年 5 月 24 日至 7 月 31 日，重复采样间隔为 28 至 66 天。

可用记录含 830 头 A2 乳型奶牛和 519 头 A1 乳型奶牛，每个牧场均包含两类。`A2_milk=1` 对应 `Betacasein=A2A2`；`A2_milk=0` 对应 `A1A1` 或 `A1A2`，与原始 `Milk_typology` 完全一致。

## 与手稿方法的关系

以奶牛作为 cluster 时，样本量为 `n=1349`，每个 cluster 有一次或两次测量。具有两次完整采样的奶牛另外组成一个分析样本。所有划分均以整头奶牛为单位。

分析的设计矩阵包含 1,060 个原始光谱变量、12 个牧场与采样批次的组合指示变量，以及 `A2_milk`，合计 1,073 列。13 个预先包含的变量已通过线性无关检查。采用两个奶牛样本组交换选择与估计角色时，每组有 674 或 675 头奶牛，设计维数超过每组奶牛数量。该数据因而直接需要高维调整。

计算使用 `ds_fit_profile` 的完整入口。每个估计样本组的变量集合由另一组执行结局搜索、profile 曲率搜索及设计回归后确定，随后在相同集合上联合计算 `h` 与 `h/2` 的系数。报告对象是包含全部光谱协变量的总体 profile 目标坐标；数值求解使用方法规定的变量集合近似这个坐标。区间由每头奶牛的联合影响值计算。

拟合前已登记 `tau=0.25,0.50,0.75`，奶牛截距采用 `Z=1` 和固定 `Lambda=1`。三个分位数的乳型坐标形成一个包含三个参数的区间族，按照手稿的 Bonferroni 规则报告同时区间。特征尺度及变量选择使用相应训练奶牛，全部目标坐标保持固定。904 头双次采样奶牛的分析报告对应总体的 profile 参数。

## 数值检查与解释条件

相邻光谱变量的相关系数中位数为 0.9986。数值输出保存选中矩阵的曲率、优化残差及目标坐标方差，最终系数在等价的加权 QR 坐标下求解。完整数值核查保存在 `results/preflight/milk_scs_20261031/`。

牧场与采样批次的组合项控制共同的采样条件。以奶牛计算区间依赖给定这些条件后的跨奶牛独立性，数据中的 6 个牧场不提供“大量独立牧场”的抽样结构。文件还缺少泌乳天数、胎次和亲缘记录，因此研究对象明确为可观测乳样特征条件下的 profile 关联。光谱在同次乳样中测量，这个条件关联的含义需要在实证章节中写清楚。

18 项完整实证拟合均已完成，每项坐标拟合时间的中位数为 17.0 秒。结果及其解释见实证结果文档。

## 复现和保存

依赖为 R 的 `readxl`、`data.table` 和 `jsonlite`，以及 Python 标准库和系统 `curl`。`readxl` 已登记到 `DESCRIPTION` 和依赖安装脚本。读取 Excel 使用 `readxl`，文件完整性检查使用标准 ZIP 库。

在仓库根目录执行以下命令。`results/scratch`、原始数据和逐条编号清单均受 `.gitignore` 管理。

```sh
mkdir -p results/scratch/runtime results/scratch/dataset_search_20261030 data/raw/candidate_datasets_20261030
curl --fail --location --output results/scratch/dataset_search_20261030/milk_metadata.json 'https://dataverse.csuc.cat/api/datasets/:persistentId?persistentId=doi:10.34810/data2338'
curl --fail --location --output data/raw/candidate_datasets_20261030/milk_README.txt https://dataverse.csuc.cat/api/access/datafile/320914
python3 research/dataset_search_20261030/download_dataset_file.py https://dataverse.csuc.cat/api/access/datafile/296867 data/raw/candidate_datasets_20261030/Database_MIR_A2_milk.xlsx 21859858
python3 research/dataset_search_20261030/verify_milk_source.py
TMPDIR="$PWD/results/scratch/runtime" Rscript research/dataset_search_20261030/audit_milk.R
TMPDIR="$PWD/results/scratch/runtime" Rscript research/dataset_search_20261030/assess_milk_design.R
```

下载命令用于缺少原始文件的环境；已有完整文件时，从 `verify_milk_source.py` 开始。汇总结果位于 `results/preflight/dataset_search_20261030/`，其中 `milk_source_manifest.json` 保存来源和版本，`milk_cohort_flow.csv` 记录样本数量，`milk_audit.json` 与 `milk_design_assessment.json` 保存核查结果。逐条编号清单位于 `data/interim/candidate_datasets_20261030/`。

仓库许可为 `CC BY-NC-ND 4.0`。GitHub 保存检索核查代码及汇总结果，原始文件使用官方获取地址。完整数据与核查文件另外备份到 `/Users/houjian/jpDMEQR-backups/dataset_search_20261030/`。
