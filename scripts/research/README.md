# Profile score 研究验证

全部命令从 `jpDMEQR` 目录运行。`TMPDIR` 使用项目内的相对路径；`JPDMEQR_CLUSTER_SOLVER=newton` 选择既有 Newton nuisance 求解函数。

```sh
export TMPDIR=results/scratch
export JPDMEQR_CLUSTER_SOLVER=newton
Rscript scripts/research/profile_score_audit_20261022.R .
Rscript scripts/research/profile_pair_audit_20261022.R .
Rscript scripts/research/profile_chart_audit_20261022.R .
Rscript scripts/research/run_profile_score_tests_20261022.R .
```

第一个入口计算重复观测的 profile score 与 Hessian、Normal 总体目标和 sandwich 极限，并在实际 cluster DGP 上完成小型优化。它还调用汇总构建脚本。第二个入口读取保存的 DGP 对象，检查配对 bandwidth 的微分关系和协方差。`profile_chart_audit_20261022.R` 使用两个独立 Normal 观察值，逐个区域计算联合密度积分，核查总概率和人口曲率，并生成 `joint_chart_table.tex`。最后一个入口运行 `test-profile-score-repeat-limit.R` 与 `test-profile-bandwidth-pair.R`。

已有完整结果时，可以使用下列命令读取结果并重建汇总。

```sh
export TMPDIR=results/scratch
export JPDMEQR_CLUSTER_SOLVER=newton
JPDMEQR_REUSE_RESEARCH_RESULTS=TRUE Rscript scripts/research/profile_score_audit_20261022.R .
JPDMEQR_REUSE_RESEARCH_RESULTS=TRUE Rscript scripts/research/profile_pair_audit_20261022.R .
Rscript scripts/research/run_profile_score_tests_20261022.R .
```

积分结果与优化对象保存在 `results/scratch/profile_score_audit_20261022/`。`repeated_geometry.csv` 保存几何检查，`population_targets.csv` 和 `population_sandwich.csv` 保存总体积分结果，`cluster_dgp.csv` 与对应 RDS 保存实际优化对象，`checks.csv` 保存检查状态。配对检查写出 `paired_components_newton.csv`、`paired_components.csv` 及 `pair_attempt_newton.csv`。测试写出 `testthat_results.rds` 和 `testthat_summary.csv`。

汇总脚本读取 `population_targets.csv`，另写 `population_target_summary.csv`。数值零使用 `numerical_zero` 标记，相应收敛阶数记录为 `NA`。

仅重建 LaTeX 表格与汇总时运行：

```sh
TMPDIR=results/scratch Rscript scripts/research/build_profile_score_summary_20261022.R .
```

表格写入工作区的 `research/profile-score-theory-20261022/population_target_table.tex`，同时写出 `population_target_table_source.csv`。表格使用指定的四个 bandwidth；显示行之间的相邻收敛阶数从误差重新计算，首行记录为 `NA`。

配对方法的完整拟合与推断检查使用下面的入口。

```sh
TMPDIR=results/scratch Rscript scripts/research/check_profile_pair_integration_20261022.R
TMPDIR=results/scratch Rscript scripts/research/run_profile_pair_registry_20261022.R .
```

注册检查采用 `config/simulation_preflight_profile_pair_20261022.csv`，运行两个设置共四次重复，核对18行方法结果与64行坐标结果。它会保存每个方法的precision交叉验证记录，并在任何方法失败、结果缺项或字段不符时终止。开发用总体目标与方向对象位于 `results/intermediate/profile_pair_20261022/`，其规模与用途由各对象的metadata记录。

当前总体对象使用 `construction_v2/` 子目录。正文与补充稿中的新表格通过下列命令从已保存的检查结果生成。

```sh
TMPDIR=results/scratch Rscript scripts/research/build_profile_pair_manuscript_20261022.R
```

默认输入为 `strict_development05` 和 `integration_20261022_05`。来源清单写入 `manuscript/generated/profile_pair_source_manifest.json`。完整计算记录见 `PROFILE_PAIR_AUDIT_20261022.md`。
