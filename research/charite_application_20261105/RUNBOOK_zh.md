# Charité 血红蛋白分析复核说明

当前研究状态保存在 `results/preflight/empirical_selection_20261101/charite_hb/study_status.json`，完整报告为同目录的 `REPORT_zh.md`。参数、患者角色和蛋白质目标已经固定，验证次数登记在 `results/preflight/empirical_selection_20261101/validation_attempts.csv`。

## 文件保存

项目目录为 `/Users/houjian/Projects/jpDMEQR`。下载文件保存在 `data/raw/covid_charite_20261105/official_cloud`，临床输入和患者角色位于 `data/derived/charite_hb_20261105`。完整模型及患者预测记录位于 `results/raw/application/charite_hb_20261105`。

三个目录分别备份到 `/Users/houjian/jpDMEQR-backups/empirical_selection_20261101` 下的 `data_raw/covid_charite_20261105`、`data_derived/charite_hb_20261105` 和 `results/application/charite_hb_20261105`。GitHub 保存研究代码、参数和汇总结果。患者记录及模型对象使用本地备份。

官方来源信息保存为 `PMC8201874.2.json` 和 `pmc8201874_listing.xml`，下载脚本读取这两个文件并核验正式版本及各文件大小。`data_manifest.json` 给出每个完整附表的官方地址。需要重新取得下载内容时，恢复来源信息后执行 `research/empirical_selection_20261101/download_charite_cloud.py`。

## 运行环境

Python 3.12.1，依赖版本列于本目录 `requirements.txt`。本次解释器为 `results/scratch/charite_env/bin/python`。R 的依赖已经列入仓库 `DESCRIPTION` 和 `scripts/00_install_dependencies.R`，具体版本记录在结果目录的 `audit_session.txt`。XGBoost 为原生3.4.1，grf 为2.6.1。

本机 XGBoost 使用已安装的 R OpenMP runtime，`DYLD_LIBRARY_PATH` 为 `/Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib`。其原生预测接口及模型保存后读取均已检查。原生库路径需要与运行机器实际安装的 OpenMP runtime 一致。

计算时将 `TMPDIR` 设为项目中的 `results/scratch/runtime`，将 `OMP_NUM_THREADS`、`OPENBLAS_NUM_THREADS` 和 `VECLIB_MAXIMUM_THREADS` 设为1，`R_MAX_VSIZE` 设为24Gb。中间文件与 Python 环境均位于 Git 忽略的 `results/scratch`。

## 从保存结果复核

在项目目录执行以下命令，可以使用现有患者预测和完整模型重新检查结果。`build_report.py` 从汇总表生成报告，表内数值来自计算结果。

```bash
results/scratch/charite_env/bin/python research/charite_application_20261105/audit_predictions.py
Rscript research/charite_application_20261105/audit_results.R
results/scratch/charite_env/bin/python research/charite_application_20261105/build_report.py
```

`audit_predictions.py` 使用独立 Python 计算重建患者损失及方法差异。`audit_results.R` 从模型系数和患者得分重建 influence、跨分位数对比及标准误。下载校验同时包括 XLSX 和 ZIP 的完整性与 PDF 解析检查。

## 计算文件顺序

`prepare_data.py` 生成患者角色和分析输入，`prepare.R` 生成训练及三个患者折的变换参数。`check_training.R` 核对求解器与训练接口，`snapshot_linear.R` 保存实现文件。线性方法通过 `research/yeast_application_20261101/run_queue.R` 执行，传入 `linear_implementation`、进程数量及当前结果目录。

`check_forest.R` 完成原生 grf 接口检查，`snapshot_forest.R` 生成两个 forest 任务，随后使用同一个队列执行程序。`export_boost.R` 导出训练折的完整精度数据。`boost.py` 在 `--mode raw` 和 `--mode offset` 下分别执行原生 XGBoost；`--check` 运行实际数据上的接口校验。

全部19项训练任务完成后，`release_validation.R` 固定蛋白质目标并登记验证。`predict_validation.R` 和 `predict_boost.py` 读取验证患者的预测变量，生成未使用验证结局的预测。`evaluate.R` 才读取结局并计算方法比较，`inference.R` 计算固定蛋白质设计下的关联对比和敏感性结果。

完整训练和推断目录具有已存在即终止的检查。保存的运行保持原样。研究复算使用保存的实现文件和同一参数；新的研究尝试需要独立目录及新的登记记录。
