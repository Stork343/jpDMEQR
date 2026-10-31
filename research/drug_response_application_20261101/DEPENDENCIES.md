# 运行依赖

数据下载使用 Python 的 `pypdl` 1.5.7 及其 `aiohttp` 依赖。HTML 文本读取使用 `beautifulsoup4`，JSON 使用 Python 标准库。临时文件保存在仓库的 `results/scratch/runtime`。

数据检查使用 R 的 `data.table` 和 `jsonlite`。作者示例使用 `fbioSelectr` 0.1.0，来源为 `https://github.com/koukoulEv/fbioSelect`，其依赖包括 `grpreg`、`doParallel`、`foreach`、`dplyr`、`iterators` 和 `splines`。作者程序保存在独立的下载目录；执行记录包含来源版本。

数据下载完成后检查文件长度，ZIP 使用标准库检查，RDS 由 R 读取。尚未完成的下载文件保留 `.part` 标识。

GDSC 原始读数的处理使用 `gdscIC50` 0.99.4，来源为 `https://github.com/CancerRxGene/gdscIC50`。数据构造另需 `matrixStats` 和 `yaml`。分位数 Lasso 比较使用项目已有的 HiGHS 程序，QRF 使用 `ranger` 0.18.0。CRAN 依赖列于 `scripts/00_install_dependencies.R`。作者程序采用 `R CMD INSTALL` 从各自下载目录安装。
