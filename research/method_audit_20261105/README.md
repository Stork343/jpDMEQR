# 临床预测方法诊断

`population_target.py` 通过确定性积分计算 profile 参数与普通分位数的差异，并检查共同常数校准后组间系数的变化。计算不生成样本。

`inspect_saved.R` 读取 RA、Pregnancy、GNHS 已保存的 90 个 JP 搜索模型，记录临床变量与蛋白质的选择情况。`fixed_support_diagnostic.R` 在原探索分组上完成四种条件的拟合，每个队列 20 次。该脚本要求新的输出目录，已存在目录会终止运行。

`verify_saved.R` 读取全部诊断拟合，复核校准最优性，并记录 ENET 的入选变量数量。`report.py` 用 Python 独立复算受试者等权损失，检查 R 结果并生成报告。

## 运行环境

工作目录为 `/Users/houjian/Projects/jpDMEQR`。使用项目已安装的 R、data.table、jsonlite、yaml、matrixStats、quantreg、Matrix、clarabel，以及既有 Rcpp 计算模块。Python 环境为 `results/scratch/charite_env/bin/python`，依赖 numpy、pandas、scipy、PyYAML 和 tabulate。R 版本及包信息见输出目录的 sessionInfo.txt。

```sh
export TMPDIR="$PWD/results/scratch/runtime"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export R_MAX_VSIZE=24Gb
results/scratch/charite_env/bin/python research/method_audit_20261105/population_target.py
Rscript research/method_audit_20261105/inspect_saved.R
Rscript research/method_audit_20261105/fixed_support_diagnostic.R config/application_ra_westchina_20261105.yml
Rscript research/method_audit_20261105/fixed_support_diagnostic.R config/application_pregnancy_detroit_20261105.yml
Rscript research/method_audit_20261105/fixed_support_diagnostic.R config/application_gnhs_baseline_bmi_20261105.yml
Rscript research/method_audit_20261105/verify_saved.R
results/scratch/charite_env/bin/python research/method_audit_20261105/report.py
```

## 输出

聚合结果和中文报告保存于 `results/preflight/method_audit_20261105`。完整模型与受试者预测保存在被 Git 忽略的 `results/raw/method_audit_20261105`。本次拟合使用源代码提交 `02a6145` 中的既有计算模块，新增诊断代码与聚合结果一同提交。

诊断固定原变量选择及 Lambda，使用各队列登记的 decision_quantile。预测保留原分位数，没有执行跨分位数排序。该设计用于理解方法组成，完整方法的独立验证仍需另行完成。
