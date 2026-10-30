# 实证计算接续

分析使用 GSE65391 的129名患者和867次就诊，分位数为0.25、0.50和0.75。预测读取已经登记的五次重复五折患者分组；每个训练集合内使用五折选择变量数量和工作惩罚。探索性基因分析读取十次患者划分及其交换角色，在独立患者中估计选择所得的三个基因。

| 计算内容 | 登记数量 |
| --- | ---: |
| 预测拟合 | 900 |
| 探索性基因分析 | 360 |
| 表达处理敏感性分析 | 12 |
| 比较方法预测拟合 | 1350 |

`inputs.R` 从当前仓库读取表达矩阵。标准处理的两个临床和模块设计已经与完成的模块计算逐项比较，最大差异为0。探针中位数处理使用每个基因对应探针的实际中位数，M5.15使用既定探针集合。训练集合决定基因筛选及标准化参数。

R1/R2数据来自NCBI登记的文件。`download_inputs.R` 调用 `pypdl` 获取文件，以文件大小和修改时间识别来源，并检查完整gzip文件。`prepare_raw.R` 使用 `limma::nec` 完成每个芯片的检测概率背景校正，offset为16，robust为FALSE。24组健康参考样本确定每个探针的R2均值与R1均值之比。`preprocessCore` 根据当前训练集合估计分位数标准化目标，随后进行log2变换。

`dispatch.R` 同时运行两个计算进程。每个预测病例和基因病例保存独立结果；表达处理敏感性分析还保存每次删除患者后的拟合。`audit_progress.R` 从保存的C3预测重新计算损失，并根据患者评分或删除患者后的系数重新计算协方差。核验表位于 `results/preflight/empirical_recovery_20261030`。

运行脚本时将 `TMPDIR` 设置为当前仓库的 `results/scratch/runtime`。R依赖通过 `scripts/00_install_dependencies.R` 安装，获取原始文件的Python依赖列于本目录的 `requirements.txt`。

```bash
Rscript research/empirical_recovery_20261030/check_inputs.R
Rscript research/empirical_recovery_20261030/check_median.R
Rscript research/empirical_recovery_20261030/start.R
Rscript research/empirical_recovery_20261030/start_downloads.R
Rscript research/empirical_recovery_20261030/complete_raw.R
Rscript research/empirical_recovery_20261030/audit_progress.R
```

调度记录保存在 `results/raw/application/empirical_recovery_20261030`。完成的结果每分钟复制到 `~/jpDMEQR-backups/empirical_recovery_20261030`。代码和仅含汇总信息的核验表提交到Git，患者层面的表达、临床和预测记录保存在忽略目录中。
