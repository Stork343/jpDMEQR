# Joint Profile Inference for High Dimensional Mixed Effects Quantile Regression

本项目保存混合效应分位数回归手稿及其研究程序，目标期刊为 Scandinavian Journal of Statistics。

当前估计方法为 `JOINT-PROFILE`。训练受试者用于选择预测变量和条件调整变量，独立的估计受试者用于拟合所选模型。两个带宽下的系数估计共同消除一阶平滑项，区间计算保留同一受试者的成对协方差。总体目标为完整 ridge-profile check loss 的最小化参数。

方法定义见 [METHOD_SPECIFICATION.md](docs/METHOD_SPECIFICATION.md)，程序入口见 [研究程序说明](research/score_validity_20261026/README_zh.md)。

## 文件位置

| 目录 | 内容 |
| --- | --- |
| `research/score_validity_20261026` | 估计方法和推断程序 |
| `research/submission_completion_20261027` | 登记实验、比较方法、结果汇总与检查 |
| `research/milk_application_20261031` | 乳样实证分析及手稿表格生成 |
| `research/joint_profile_manuscript_20261026` | 固定设计计算和生成程序 |
| `R` | 数据生成、数值计算及比较方法接口 |
| `config` | 实验设计、随机种子和数值设置 |
| `manuscript` | 手稿和 Supporting Information |
| `data` | 研究数据及其来源记录 |
| `results` | 逐次计算结果、检查记录和汇总材料 |

本地工作目录为 `/Users/houjian/Projects/jpDMEQR`。GitHub 计算分支为 `recovery/research-20261028`。

## 环境准备

在项目目录执行：

```bash
mkdir -p results/scratch/runtime
export TMPDIR="$PWD/results/scratch/runtime"
Rscript scripts/00_install_dependencies.R
Rscript research/score_validity_20261026/setup.R
```

比较程序使用项目目录中的 HiGHS 源代码安装版本，最低版本为1.15.1-1。依赖安装程序同时检查 fastclime 1.4.1.1。

## 数值检查

```bash
Rscript research/score_validity_20261026/validate_joint.R
Rscript research/score_validity_20261026/validate_weighted_profile.R
Rscript research/submission_completion_20261027/check_profile_charts.R
Rscript research/joint_profile_manuscript_20261026/validate_fixed.R
```

总体目标文件按 `config/joint_profile_submission_20261027.yml` 定位。正式使用前，运行 `research/submission_completion_20261027/validate_assets.R`，检查各文件的实验身份、独立重复数量和数值精度。

## 登记计算

实验设置和重复次数由 `config/simulation_main.csv` 给出。方法映射与随机种子来自 `config/joint_profile_submission_20261027.yml`，比较方法的数值设置保存在 `config/submission_comparator_numerics_20261028.yml`。

准备好总体目标后，执行预实验并读取其检查结果：

```bash
Rscript research/submission_completion_20261027/launch_simulation.R pilot joint 8
Rscript research/submission_completion_20261027/summarise_pilot.R
```

预实验和总体目标检查通过后运行正式计算：

```bash
Rscript research/submission_completion_20261027/launch_simulation.R main joint 8
```

已有计算目录使用 `resume_simulation.R` 继续处理尚未完成的重复。完整队列通过 `collect_simulation.R` 读取全部登记编号，核对方法记录后生成汇总输入。

## 保存与备份

每次重复完成后，通过临时文件写入并改名的方式保存完整 RDS。计算进程每完成100次重复，以及完成批次时，将完整结果复制到 `~/jpDMEQR-backups`。`JPDMEQR_BACKUP_DIR` 可以指定备份位置。

`research/submission_completion_20261027/backup_live_results.R` 每60秒检查进度并复制新增的完整结果，直至全部登记队列完成。

源码和手稿通过 Git 提交并推送。研究数据与大型计算对象保存在本地，源码中的结果生成程序负责构建表格和图。运行与保存说明见 [LOCAL_EXECUTION_20261028.md](docs/LOCAL_EXECUTION_20261028.md)。

## 手稿与数据

手稿入口为 `manuscript/sjs_profile_quantile_main.tex`，Supporting Information 入口为 `manuscript/sjs_profile_quantile_supplement.tex`。图表生成程序读取已经保存的计算结果。

乳样实证已经完成 18 项登记计算。数据来源见 [乳样核查](research/dataset_search_20261030/README_zh.md)，完整分析见 [实证方案](research/milk_application_20261031/PROTOCOL_zh.md)和 [实证结果](research/milk_application_20261031/RESULTS_zh.md)。比较方法的原始参考程序保存在 `data/raw/public_methods`，其来源说明位于 `docs/benchmark_notes`。

软件许可为 GPL version 3 或更高版本，引用信息见 `inst/CITATION`。
