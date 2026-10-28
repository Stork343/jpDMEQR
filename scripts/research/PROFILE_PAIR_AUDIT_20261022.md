# Profile pair 计算核查

## 几何计算

`results/scratch/profile_pair_geometry_development_20261022/strict_development05/manifest.json` 记录20个 deterministic seeds，覆盖四个 tau、q=1/2、三个 h 倍率及不等 cluster sizes。81项 geometry 对象来自 `strict_development04`，其中20项使用 non-diagonal Lambda。05另行完成20项真实 fit/infer，每项请求规范中的八个 Dantzig 倍率，并保存全部 CV 对象。

gradient 最大误差为 `2.6907769e-8`，Hessian 最大误差为 `2.18181298e-4`；对应接受阈值为 `2e-4` 与 `3e-3`。profile identity、sandwich、cross-bandwidth covariance 与 assembly 的最大误差分别为 `6.66134e-16`、`5.20417e-18`、`2.77556e-16` 和 `2.66454e-15`。nuisance gradient 最大值为 `9.93597e-8`。

79项几何用例的 Dantzig 求解可行。两项不可行用例通过独立 ECOS minimax-defect 计算核查，其最小可行误差分别为 `0.0848864713` 与约 `1`，大于固定诊断阈值 `0.05`。负曲率用例的 paired Hessian 特征值为 `-0.3733592`，微分及同 score sandwich 关系通过。

05的运行设置为 `clarabel_native/newton/sequential/CV1/Hrefresh1/sharedFALSE`。04的环境字段使用 `native`；两个入口显式调用同一 `CLARABEL-native` 函数。04保存的环境记录保持原值。

## 总体对象

`R/v2_simulation_v2_part04.R` 提供 `profile-pair-population-assets-20261022-v2`，保存位置为 `results/intermediate/profile_pair_20261022/construction_v2/{targets,directions}/<experiment_id>.rds`。验证函数要求当前 score_definition 与 construction_version，final 模式还要求真实 implementation commit、`n_population>=100000`、`repeats>=4` 及通过数值接受条件。

target 对象保存独立重复的 paired roots、target/refinement/analysis bandwidth、MCSE 和数值检查。direction 在 structural 设置使用生成器的 beta0，并检查每次生成的 beta0 完全相同；profile_mc 设置使用独立 target 对象的 beta_star_mc。每个对象及每次重复均保存 target_source 和 evaluation_point。

实际 development 构建标识为 `PAIR-ASSET-DEV-structural_direction_check01`，使用120个 clusters、两次独立重复，已完成写入与回读。direction evaluation_point 为 `(0.75,0)`，独立 target MC root 为 `(0.67390736,-0.01436219)`。target 最大 score 残差为 `1.718134e-8`，最大 nuisance gradient 为 `9.954257e-8`，covariance identity 误差为 `4.857226e-17`，inverse defect 为 `1.781369e-17`。该对象的 `development=TRUE`，`final_eligible=FALSE`。

## 结果记录与回归检查

failed paired inference 为每个预设坐标保存 NA 估计、NA 区间与 failed 状态。覆盖率汇总同时提供包含失败的分母和仅含可用区间的统计。`p_value` 使用零假设 beta=0；`wald_z` 使用研究目标作为中心。paired Bahadur 与 L_k 使用完全匹配 beta、h 和 design 的总体 score 均值，匹配条件不满足时保存 NA 与 reason。

真实保存结果 `integration_20261022_02` 中的三项失败重建出12条缺少的坐标记录，完整记录共64条。重建与汇总文件保存在 `results/scratch/profile_pair_failure_accounting_20261022/`，原始结果文件保持原值。

`test-profile-pair-geometry-audit.R` 有18项通过的断言，`test-profile-pair-population-assets.R` 有40项，`test-profile-pair-runner-accounting.R` 有30项。后者使用实际 DGP 检查失败记录、覆盖率分母、p_value、总体 centering、Bahadur、L_k 和全部 NA 的汇总。断言对象分别保存在 `results/scratch/profile_pair_geometry_tests_20261022/`、`results/scratch/profile_pair_asset_tests_20261022.rds` 与 `results/scratch/profile_pair_accounting_tests_20261022/`。

全部入口从 repository 根目录运行，使用相对路径 `TMPDIR=results/scratch`。

```sh
TMPDIR=results/scratch Rscript scripts/research/audit_profile_pair_geometry_20261022.R . <new_run_id> strict_development04
TMPDIR=results/scratch Rscript scripts/research/build_profile_pair_assets_20261022.R . <new_run_id> structural
TMPDIR=results/scratch Rscript scripts/research/run_profile_pair_geometry_tests_20261022.R
TMPDIR=results/scratch Rscript scripts/research/run_profile_pair_accounting_tests_20261022.R
```
