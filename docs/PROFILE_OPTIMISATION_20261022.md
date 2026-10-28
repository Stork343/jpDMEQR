# Profile数值计算

本轮计算使用`compiled` cluster nuisance引擎。完整criterion、ridge矩阵Lambda及原始coordinate的KKT和coefficient变化条件由`METHOD_SPECIFICATION.md`规定。数值配置通过`solver_configuration_v2()`进入运行记录。

## Unpenalised Newton

`JPDMEQR_UNPENALISED_OPTIMIZER=newton_armijo`仅在全部`lambda_beta * penalty_factor`等于零时启用。每次函数计算都重新求解完整cluster nuisance问题。内部stationarity tolerance为`min(requested_grad_tol, 1e-11)`，两项数值分别保存。

对coordinate j，固定计算尺度为

\[
s_j^2=\frac1n\sum_{i=1}^n\frac1{m_i}\sum_{l=1}^{m_i}X_{ilj}^2,
\qquad D=\operatorname{diag}(s_1,\ldots,s_p).
\]

Newton方程使用`D^{-1} H D^{-1}`与`D^{-1} g`。其特征值在求解方程时设下限`1e-8 * max(1, max(abs(eigenvalues)))`，方向转换回原始coefficient单位。一次方向上的最大预测量变化限制在`max(1, max(abs(y-X beta)))`以内。Armijo条件使用完整profile objective，常数为`1e-4`，步长每次乘以`0.5`。浮点数计算容差为`64 * machine_epsilon * max(1, abs(objective))`。

每次更新后检验原始score。接受条件仍包含指定KKT tolerance和原始coefficient的相对变化条件。最终重新计算完整Hessian、两种Hessian表达式及cluster score，失败继续终止当前调用。未指定该数值配置时使用`proximal_gradient`。

## Cluster RMS proximal

`JPDMEQR_PROXIMAL_SCALING=cluster_rms`使用同一固定D。令`theta=D beta`，转换后的光滑部分Hessian为`D^{-1} H D^{-1}`，L1 penalty为`lambda sum_j penalty_factor_j abs(theta_j)/s_j`。转换回原始单位的更新为

\[
\beta_j^{new}=S\left(\beta_j-tg_j/s_j^2,
                      t\lambda\,penalty\_factor_j/s_j^2\right).
\]

line search的二次项为`sum((delta * s)^2)/(2t)`。完整目标函数中的penalty与原始coordinate KKT保持原定义；momentum及相对coefficient变化条件使用原始beta。`none`使用单位尺度。

## 数值检查

`scripts/research/check_profile_unpenalised_newton_20261022.R`对真实129名subject、867次就诊的A16与B17设计计算三个既定tau和q1/q2。每个Newton解与原proximal函数在相同criterion下比较原始coefficient、objective和KKT。机器可读比较结果保存于版本化的`results/scratch/unpenalised_newton_20261022_v3/`。

production配置的DGP检查覆盖q1/q2、非对角Lambda和三个tau。`test-profile-proximal-scaling.R`另外核查mixed penalty的共同解与原始KKT。真实prediction训练检查保存完整5000gene筛选及100gene设计，在固定最终penalty下比较两条数值路径。

## 完整workflow存储

`checkpoint_storage=compact`保留每次preprocessing参数、subject与sample列表、probe映射及module audit。clinical formula只引用数据中的列名，其环境设为`baseenv()`，计算矩阵由共同输入和这些参数重建。各次coefficient、迭代history、KKT、nuisance诊断及Hessian曲率诊断独立保存。完整样本的cluster score和Hessian，以及全部leave-one coefficient矩阵保存在workflow artifact中。

`read_compact_application_workflow_v2()`读取全部checkpoint，重新组成jackknife矩阵与两种区间，要求与保存的完整表格精确相同。matched runner进一步转换全部coefficient至固定总体单位，核查重读后的完整报告表格。共同输入由固定版本路径、构建commit、生成seed和抽样subject身份记录引用。
