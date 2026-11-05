import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

root = Path("results/preflight/method_audit_20261105")
configs = [yaml.safe_load(path.read_text()) for path in [
    Path("config/application_ra_westchina_20261105.yml"),
    Path("config/application_pregnancy_detroit_20261105.yml"),
    Path("config/application_gnhs_baseline_bmi_20261105.yml"),
]]
labels = {"ra_westchina_20261105": "RA", "pregnancy_detroit_20261105": "Pregnancy",
          "gnhs_baseline_bmi_20261105": "GNHS"}
variants = ["profile_selected", "pinball_selected", "profile_full_clinical", "pinball_full_clinical"]
all_checks = []
loss_records = []
maximum_loss_error = 0.0
for config in configs:
    study = config["study_id"]
    frame = pd.read_csv(Path("results/raw/method_audit_20261105") / study / "predictions.csv")
    assert not frame.duplicated(["sample_id", "variant"]).any()
    assert frame.groupby("strain_id")["fold"].nunique().eq(1).all()
    assert frame.groupby("variant")["sample_id"].nunique().nunique() == 1
    residual = frame["y"] - frame["prediction"]
    loss = residual * (frame["tau"] - (residual < 0))
    maximum_loss_error = max(maximum_loss_error, np.abs(loss - frame["loss"]).max())
    frame["independent_loss"] = loss
    patient = frame.groupby(["variant", "strain_id"])["independent_loss"].mean()
    calculated = patient.groupby("variant").mean()
    saved = pd.read_csv(root / study / "fixed_support_losses.csv").set_index("variant")
    maximum_loss_error = max(maximum_loss_error, np.abs(calculated - saved["loss"]).max())
    baseline = calculated["profile_selected"]
    record = {"队列": labels[study], "受试者": frame["strain_id"].nunique(),
              "tau": config["decision_quantile"], "原目标原变量": baseline}
    changes = pd.read_csv(root / study / "fixed_support_comparisons.csv").set_index("variant")
    for variant, label in zip(variants[1:], ["取消患者效应", "保留全部临床变量", "同时调整"]):
        relative_change = calculated[variant] / baseline - 1
        assert abs(relative_change + changes.loc[variant, "relative_gain"]) < 1e-10
        record[label] = relative_change
    loss_records.append(record)
    checks = pd.read_csv(root / study / "solver_checks.csv")
    assert len(checks) == 20 and checks["gap"].abs().max() < 1e-8
    assert checks["beta_kkt"].max() < 1e-8 and checks["gamma_kkt"].max() < 1e-8
    assert checks["original_prediction_error"].dropna().max() < 1e-6
    all_checks.append(checks)
assert maximum_loss_error < 1e-10
solver = pd.concat(all_checks, ignore_index=True)
calibration = pd.read_csv(root / "calibration_checks.csv")
assert len(calibration) == 60
assert calibration["solver_status"].isin(["Solved", "AlmostSolved"]).all()
assert (calibration["below"] <= calibration["tau"] + 1e-10).all()
assert (calibration["at_or_below"] >= calibration["tau"] - 1e-10).all()
warnings = pd.read_csv(root / "calibration_warnings.csv")
assert warnings["message"].eq("Solution may be nonunique").all()
selection = pd.read_csv(root / "saved_selection_audit.csv")
assert len(selection) == 90
enet = pd.read_csv(root / "enet_selection.csv")
population = pd.read_csv(root / "population_targets.csv")
slopes = pd.read_csv(root / "population_slopes.csv").set_index("case")
probability = pd.read_csv(root / "population_calibration.csv")
population_check = json.loads((root / "population_checks.json").read_text())
assert population_check["status"] == "passed" and population_check["random_samples_generated"] == 0
normal = population.query("tau == 0.75 and Lambda == 1").copy()
normal = normal[["cluster_size", "normal_quantile", "profile_target", "achieved_probability"]]
normal.columns = ["每人记录数", "分位数", "profile 参数", "对应累计概率"]
normal["每人记录数"] = normal["每人记录数"].astype(int).astype(str)
for column in ["分位数", "profile 参数", "对应累计概率"]:
    normal[column] = normal[column].map(lambda value: f"{value:.6f}")
diagnostic = pd.DataFrame(loss_records)
diagnostic.to_csv(root / "diagnostic_summary.csv", index=False)
display = diagnostic.copy()
for column in ["取消患者效应", "保留全部临床变量", "同时调整"]:
    display[column] = display[column].map(lambda value: f"{value:+.2%}")
display["原目标原变量"] = display["原目标原变量"].map(lambda value: f"{value:.6f}")
display["tau"] = display["tau"].map(lambda value: f"{value:.2f}")
decision = selection[selection["decision_quantile"] & selection["method"].eq("JOINT_DIRECT")]
capacity = []
for study, group in decision.groupby("study_id"):
    other = enet[enet["study_id"].eq(study)]
    capacity.append({"队列": labels[study], "JP 入选蛋白质范围":
                     f"{group['proteins_selected'].min()} 至 {group['proteins_selected'].max()}",
                     "达到搜索上限的折数": int((group["selected_additions"] == group["maximum_additions"]).sum()),
                     "ENET 入选蛋白质范围":
                     f"{other['selected_proteins'].min()} 至 {other['selected_proteins'].max()}"})
capacity = pd.DataFrame(capacity)
checks = {"status": "passed", "saved_models_checked": len(solver), "saved_search_models_audited": len(selection),
          "maximum_independent_loss_error": float(maximum_loss_error),
          "maximum_original_prediction_error": float(solver["original_prediction_error"].max()),
          "maximum_primal_dual_gap": float(solver["gap"].abs().max()),
          "maximum_beta_kkt": float(solver["beta_kkt"].max()),
          "maximum_gamma_kkt": float(solver["gamma_kkt"].max()),
          "calibration_nonunique_warnings": len(warnings),
          "calibration_optimality_checks_passed": len(calibration),
          "formal_validation_used": False, "manuscript_adoption": False}
(root / "verification.json").write_text(json.dumps(checks, indent=2) + "\n")
shape = slopes.loc["conditional_shape_repeated"]
shape_probability = probability[probability["case"].eq("conditional_shape_repeated")].sort_values("group")
report = rf"""# 方法目标与临床预测诊断

## 估计目标与预测评价

现行方法估计固定 Lambda 下的 regularised profile 参数。新受试者的预测使用固定效应与一个训练样本校准常数，评价标准为受试者等权 pinball loss。这两个目标通常具有不同的总体最优解。该差异在样本数增加后仍可存在，也可能改变协变量的系数。

方法规范定义

\[
\beta^*(\Lambda)=\arg\min_\beta E\min_\gamma\left\{{
\frac1m\sum_j\rho_\tau(Y_j-X_j^T\beta-Z_j^T\gamma)
+\frac12\gamma^T\Lambda\gamma\right\}}.
\]

当每人一条记录且 Z=1 时，写 r=Y-X^T beta，患者效应为

\[
\widehat\gamma=\operatorname{{clip}}\left(r,\frac{{\tau-1}}{{\lambda}},\frac{{\tau}}{{\lambda}}\right),
\qquad \psi_{{prof}}(r)=\operatorname{{clip}}(\lambda r,\tau-1,\tau).
\]

因此总体方程包含经过截断的残差。普通 quantile regression 的总体方程包含残差是否小于零的指示函数。下面对标准正态分布直接积分，tau=0.75，lambda=1。重复记录在给定设计后独立，这是现行模型允许的情形。总体计算省略随样本数增加而趋于零的预测 ridge，lambda 保持固定。

{normal.to_markdown(index=False, disable_numparse=True)}

单个分布的截距差可以通过总体校准消除。协变量对应的分布形状改变时，一个共同常数无法修复所有条件分位数。考虑两组受试者，每人两条记录，两组等权。第一组的分布为 0.8 N(0,0.05²)+0.2 N(4,0.05²)，第二组为 N(0,1)。固定设计包含截距和组别，Z=1，lambda=1，tau=0.75。真实分位数组间差为 {shape['quantile_slope']:.6f}，profile 系数为 {shape['profile_slope']:.6f}。加入最优共同校准常数后，组间差仍为 {shape['calibrated_slope']:.6f}；两组累计概率分别为 {shape_probability.iloc[0]['conditional_probability']:.2%} 和 {shape_probability.iloc[1]['conditional_probability']:.2%}，平均恰为 75%。这些数值来自确定性积分，没有生成临床样本。

对于 m 条独立记录，令 c_j=(tau-1+j/m)/lambda，患者效应是 m 个残差与 m+1 个 c_j 的中位数。在 c_j 与 c_(j+1) 之间，其分布函数为 P[Binomial(m,F_Y(a+t)) ≥ m-j]。脚本对该分布积分求 E gamma=0，并用单记录闭式得分及总体目标的数值导数交叉检查。

手稿的 paired coefficient 构造消除 bandwidth 展开的一阶项，目标仍是 beta*(Lambda)。当前预测比较使用 unsmoothed JOINT_DIRECT 和 JOINT_OFFSET 搜索模型。预测损失可以评价这两个分支的应用价值；paired coefficient 的区间覆盖与参数估计效率需要各自对应的检验。

## 临床变量的处理

`docs/METHOD_SPECIFICATION.md` 的 Prediction for new subjects 要求保留全部预定临床与 module 变量。`research/clinical_search_20261105/panel.R` 中 cp_features 的 penalty_factor 只固定截距，其余列都参与筛选。保存的 90 个 JP 搜索模型已逐一读取。

在登记的 decision_quantile，RA 和 GNHS 的 JOINT_DIRECT 及 JOINT_OFFSET 各五折都省略了部分临床列。Pregnancy 的设计只有蛋白质与截距，这项差异不涉及该数据集。JOINT_OFFSET 的 ridge 背景仍使用临床信息；这里统计的是其分位数修正部分能否自由估计完整临床系数。

搜索所得变量数量也值得检验。

{capacity.to_markdown(index=False)}

这些计数来自同一探索分组的保存模型。ENET 的数量按绝对系数大于 1e-9 计数。模型容量差异已经确认，其对预测差距的贡献需要重新训练并比较，单靠计数无法确定。

## 固定变量集合后的诊断重算

三个队列沿用既有五折受试者分组，每个队列计算四种条件，共 60 次拟合。比较固定原入选变量时保留或取消患者效应，并检查加入完整临床变量后的对应结果。所有条件保留原训练预处理、预测 ridge 和训练校准。原方法重算预测与保存结果的最大差异为 {checks['maximum_original_prediction_error']:.2e}。

下表采用各队列登记的 decision_quantile。百分比为相对原方法的损失变化，负数表示改善。每位受试者的记录内平均，再对受试者等权平均。比较保持该分位数的原始预测，不执行跨分位数排序；因此绝对损失与以排序预测生成的既有筛选表分别报告。

{display.to_markdown(index=False, disable_numparse=True)}

取消患者效应没有改善这三个队列的固定变量预测。补入全部临床变量后，GNHS 的改善幅度有限，RA 的损失增加。Pregnancy 没有额外临床变量，两组对应计算相同。该诊断没有重新选择变量或 Lambda；结果描述的是在既有选择条件下改变拟合目标和临床设计的作用。

所有保存结果均通过独立 Python 损失复算，最大差异为 {checks['maximum_independent_loss_error']:.2e}。60 次拟合的最大 primal dual gap 为 {checks['maximum_primal_dual_gap']:.2e}。校准重算出现 {checks['calibration_nonunique_warnings']} 次 Solution may be nonunique 提示；60 次校准均满足加权分位数最优性条件。详细记录见 calibration_warnings.csv 和 calibration_checks.csv。

## 研究设计需要回答的问题

| 问题 | 已完成的检查 | 当前判断 | 后续检验 |
|:--|:--|:--|:--|
| profile 参数是否就是临床条件分位数系数 | 总体积分与解析得分检查 | 一般情形下不同，共同校准也可能保留系数差异 | 明确临床参数，再证明所用估计方程与该参数一致 |
| 未保留临床变量能否解释预测不足 | 保存模型读取与完整临床设计重算 | 实现与规范存在差异，补入变量没有普遍改善 | 在完整训练选择流程中固定预定临床设计 |
| 固定变量后取消患者效应能否改善预测 | 三个队列共 60 次拟合 | 三个队列均未改善 | 在同一科学目标下比较利用记录相关性的方法 |
| 搜索容量是否限制预测 | JP 与 ENET 保存系数检查 | 容量不同，影响尚未定量确定 | 对共同输入和训练分组比较完整选择程序 |
| paired coefficient 推断是否有效 | 本次未重跑区间实验 | 本次预测诊断没有评价其覆盖率 | 按对应参数报告覆盖率与区间长度 |

现阶段重设计应围绕明确的临床参数展开。若应用要求预测新受试者的条件分位数，应把受试者等权条件分位数风险写入估计目标，并研究如何利用同一受试者内的相关性改善估计效率。患者效应的处理需要给出保持该参数的推导；Lambda 不能同时承担决定科学参数和改善预测的两种职责而缺少连接条件。

新的完整比较需要保留预定临床变量，统一受试者分组，并检验不同变量选择程序对稀疏与分散信号的适应能力。参数由训练队列内的受试者交叉验证决定。修改后的程序需要在未参与方法开发的受试者上评价，才能判断预测优势能否重复出现。

本次结果属于探索诊断。三个队列保留的验证受试者没有用于本次计算，formal_validation_used=false，manuscript_adoption=false。

## 计算文件与来源

数学定义来自 `docs/METHOD_SPECIFICATION.md`。预测流程来自 `research/empirical_selection_20261101/prediction_core.R`，特征接口为 `research/clinical_search_20261105/panel.R`。原始模型路径由三个 application 配置中的 output_root 指定，诊断拟合保存在 `results/raw/method_audit_20261105`。

数据来源及清理过程见 `results/preflight/empirical_selection_20261101/ra_westchina/REPORT_zh.md`、`results/preflight/empirical_selection_20261101/pregnancy_detroit/REPORT_zh.md` 和 `results/preflight/empirical_selection_20261101/gnhs_baseline_bmi/REPORT_zh.md`。本次只使用已有探索受试者及其保存模型。

`research/method_audit_20261105` 保存总体积分、模型读取、诊断拟合与复核代码。表格由这些程序输出生成，verification.json 保存全部数值检查。
"""
(root / "REPORT_zh.md").write_text(report.strip() + "\n")
print(json.dumps(checks, indent=2))
