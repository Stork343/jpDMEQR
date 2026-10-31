# 运行依赖

R 使用项目已有的 `data.table`、`jsonlite`、`yaml`、`glmnet`、`quantreg`、`RcppEigen`、`clarabel` 和 `processx`。分位数森林调用 `grf` 2.6.1；安装入口为 `scripts/00_install_dependencies.R`。精确分位数 Lasso 使用项目的 HiGHS 1.15.1 编译程序。

Python 的 `pandas` 和 `beautifulsoup4` 用于数据检查及来源读取。IBDMDB 的 BIOM 文件使用 `biom-format` 读取。下载程序使用 `pypdl`。环境与安装缓存保存在 `results/scratch`，运行的 `TMPDIR` 指向 `results/scratch/runtime`。

源文件和配置提交至 GitHub，原始数据及逐次拟合对象另行保存于 `/Users/houjian/jpDMEQR-backups/empirical_selection_20261101`。
