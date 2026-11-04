# 分析环境

R 使用 BEDMatrix 2.0.4、data.table、jsonlite、yaml、glmnet 4.1-9、quantreg、clarabel 及项目已经验证的 HiGHS quantile Lasso 求解器。PLINK 2.0.0-b.1.3 根据训练菌株的 SNP 频率与 LD 构造输入。

Python 环境为 `results/scratch/efa_env`，使用 pandas、numpy、PyYAML、xlrd 2.0.2 和 bed_reader 1.1.0。xlrd 读取原作者发布的 Excel 文件，bed_reader 与 R 的 BEDMatrix 交叉检查真实基因型。

执行时将 TMPDIR 设置为项目的 `results/scratch/runtime`，计算过程限定 BLAS 和 OpenMP 的线程数为 1。大矩阵分析使用 R_MAX_VSIZE=24Gb。程序在输入、计算或最优性检查失败时终止。
