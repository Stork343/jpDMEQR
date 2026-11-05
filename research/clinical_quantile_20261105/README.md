# 临床条件分位数计算

当前点估计使用 CENTERED_SMOOTHING.md 中的中心化平滑，CQ ENET 提供非平滑比较。METHOD.md 定义患者正交得分和协方差加权组合。THEORY_REVIEW.md 记录已经核验的推导及收敛速度要求，RELATED_WORK.md 保存相近研究的比较。

当前运行配置为 `config/clinical_quantile_smooth_20261106.yml`。启动 `run_smooth_queue.py` 后，三个进程按患者折计算，并保存每个训练参数路径。全部队列的验证受试者继续保留。

```sh
mkdir -p results/scratch/runtime
export TMPDIR="$PWD/results/scratch/runtime"
results/scratch/charite_env/bin/python -m venv results/scratch/cq_solver_env
results/scratch/cq_solver_env/bin/python -m pip install -r research/clinical_quantile_20261105/requirements-solver.txt
results/scratch/charite_env/bin/python research/clinical_quantile_20261105/run_smooth_queue.py
results/scratch/charite_env/bin/python research/clinical_quantile_20261105/summarise.py config/clinical_quantile_smooth_20261106.yml
Rscript research/clinical_quantile_20261105/verify_predictions.R config/clinical_quantile_smooth_20261106.yml
```

## 比较方法运行

工作目录为 `/Users/houjian/Projects/jpDMEQR`。R 依赖均已列入项目安装脚本，包括 quantreg、clarabel、osqp、glmnet、Rcpp、Matrix、data.table、yaml、jsonlite 和 matrixStats。HiGHS 使用已有的 `load_pooled_lp.R` 接口。Python 使用 `results/scratch/charite_env/bin/python`，需要 numpy、pandas、scipy、PyYAML 和 tabulate。

```sh
export TMPDIR="$PWD/results/scratch/runtime"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export R_MAX_VSIZE=24Gb
Rscript research/clinical_quantile_20261105/check_engine.R
Rscript research/clinical_quantile_20261105/check_orthogonal.R
Rscript research/clinical_quantile_20261105/check_inference.R
```

`run_study.R` 接收一个既有 application 配置，使用其中已保存的患者分组和训练预处理。输出目录存在时会终止，不覆盖既有运行。

```sh
Rscript research/clinical_quantile_20261105/run_study.R config/application_ra_westchina_20261105.yml
Rscript research/clinical_quantile_20261105/run_study.R config/application_pregnancy_detroit_20261105.yml
Rscript research/clinical_quantile_20261105/run_study.R config/application_gnhs_baseline_bmi_20261105.yml
```

完成一个队列后，可以汇总并独立复算该队列。未完成队列不会进入结果比较。

```sh
results/scratch/charite_env/bin/python research/clinical_quantile_20261105/status.py
results/scratch/charite_env/bin/python research/clinical_quantile_20261105/summarise.py
Rscript research/clinical_quantile_20261105/verify_predictions.R
Rscript research/clinical_quantile_20261105/run_inference.R config/application_ra_westchina_20261105.yml C_age
results/scratch/charite_env/bin/python research/clinical_quantile_20261105/aggregate_inference.py
Rscript research/clinical_quantile_20261105/verify_inference.R
results/scratch/charite_env/bin/python research/clinical_quantile_20261105/write_report.py
```

GNHS 的点估计完成后，`run_inference.R` 不带参数执行配置中的 C_prior_bmi 检查。预先确定的 clinical 列用于计算方法的探索检验，正式应用的科学问题需要单独登记。

## 保存

配置为 `config/clinical_quantile_20261105.yml`。完整拟合在 `results/raw/clinical_quantile_qp_20261105`，聚合结果在 `results/preflight/clinical_quantile_20261105`。原始患者记录和模型对象不进入 Git。已完成的队列与对应源代码另外保存到 `/Users/houjian/jpDMEQR-backups/empirical_selection_20261101`。

当前采用条件分位数目标。原 profile 方法定义及其模拟保持独立的版本记录。各预测队列的验证受试者继续保留。
