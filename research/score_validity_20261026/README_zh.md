# 联合剖面估计与区间推断

当前方法入口为 `ds_fit_profile`，标识为 `joint_profile_design_20261026`。它使用条件相关性筛选需要调整的协变量，并在独立样本组上联合估计两个带宽的完整profile模型。估计对象仍是固定 `Lambda` 下的总体profile系数。

[独立验证结果](RESULTS.md)记录了1,440次模拟和计算时间。[方法定义](METHOD.md)给出全部计算规则，相关误差界与目标参数结论见[理论推导](JOINT_THEORY.md)。[文献比较](NOVELTY.md)说明算法贡献与已有工作的联系。

## 运行示例

工作目录使用仓库根目录 `jpDMEQR`。R环境需要 `RcppEigen`、`glmnet` 和 `filelock`，依赖安装脚本为 `scripts/00_install_dependencies.R`。`setup.R`把当前仓库的包安装到本方法的独立运行目录。

```sh
mkdir -p results/scratch
TMPDIR=results/scratch Rscript research/score_validity_20261026/setup.R
TMPDIR=results/scratch Rscript research/score_validity_20261026/example_joint.R
```

在R中使用自己的数据时，调用方式如下。`cluster_id`标记完整的独立样本组，`X`需要唯一的列名。目标坐标可以使用列名或整数位置。

```r
source("research/score_validity_20261026/load_joint.R")
fit <- ds_fit_profile(y, X, Z, cluster_id,
                      tau = 0.5,
                      coordinates = c("x00001", "x00002"),
                      Lambda = diag(ncol(Z)))
fit$table
```

`fit$table`包含估计值、标准误和置信区间。两种带宽的联合拟合及其cluster贡献保存在 `fit$details`，变量集合与样本角色保存在 `fit$training`。默认置信水平为95%。

## 复现模拟

以下命令使用四个R进程运行 `evaluation.csv` 的24种设置，每种60次。输出目录必须使用新的名称，已有原始结果会受到保护。

```sh
TMPDIR=results/scratch Rscript research/score_validity_20261026/run_validation.R reproduction_run
```

随机种子按设置和重复编号确定。每个任务在计算前保存请求，发生错误时保留失败记录并终止执行。编译缓存使用进程锁。

已有独立验证的原始目录是 `results/scratch/score_validity_20261026/design_independent`。覆盖率与逐坐标统计位于 `results/tables/score_validity_20261026/design_independent`。下列命令重新生成汇总和目标参数检查。

```sh
TMPDIR=results/scratch Rscript research/score_validity_20261026/summarise_evaluation.R design_independent
TMPDIR=results/scratch Rscript research/score_validity_20261026/audit_selected_targets.R design_independent
TMPDIR=results/scratch Rscript research/score_validity_20261026/build_report.R design_independent
```

## 计算规则

变量搜索检查全部 `p` 个协变量。最终系数拟合只求解已选变量的小矩阵，每个cluster的nuisance方程维数为 `q`。时间报告对应四个目标系数，增加目标坐标会增加相应的条件回归和拟合计算。

`penalty_factor`使用零值指定必须包含的协变量，一值表示允许搜索的协变量。数值计算采用完整profile梯度及Schur曲率，并检查优化残差。数据生成参数只出现在模拟与结果检查脚本中。

## 文件入口

`small_profile.R`实现选中模型的凸优化，`joint_profile_inference.R`构造最终估计和cluster方差。`method.R`提供变量搜索与曲率筛选函数，`geometry.cpp`计算所需的小矩阵。数值检查使用 `validate_weighted_profile.R`，估计接口检查使用 `validate_joint.R`。
