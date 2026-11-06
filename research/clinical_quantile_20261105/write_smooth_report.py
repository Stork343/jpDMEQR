from pathlib import Path
import json

import numpy as np
import pandas as pd

base = Path("results/preflight")
out = base / "clinical_quantile_smooth_20261106"
current = base / "clinical_followup_20261106"
current.mkdir(exist_ok=True)
cohorts = {"RA": ("ra_westchina_20261105",.75), "Pregnancy": ("pregnancy_detroit_20261105",.5),
           "GNHS": ("gnhs_baseline_bmi_20261105",.75)}
decisions, intervals = [], {}
for label, (study, tau) in cohorts.items():
    folder = out if label == "Pregnancy" else base / f"cq_{study}_pw1"
    x = pd.read_csv(folder/"discovery_comparisons.csv")
    x = x.loc[x.study_id.eq(study) & x.tau.eq(tau)].set_index("comparator")
    assert len(x) >= 2 and x.index.is_unique
    np.testing.assert_allclose(x.relative_improvement, 1-x.current_loss/x.comparator_loss, rtol=0, atol=1e-12)
    assert json.loads((folder/"prediction_acceptance.json").read_text())["status"] == "passed"
    best = x.loc[x.comparator_loss.idxmin()]
    decisions.append({"study_id": study, "tau": tau, "best_comparator": best.name, "new_loss": best.current_loss,
                      "best_comparator_loss": best.comparator_loss, "relative_gain_over_best": best.relative_improvement,
                      "relative_gain_over_joint_direct": x.at["JOINT_DIRECT","relative_improvement"] if "JOINT_DIRECT" in x.index else None,
                      "linear_screen_passed": bool(x.relative_improvement.min()>=.02), "formal_validation_used": False,
                      "manuscript_adoption": False, "compared_methods": "|".join(x.index), "authoritative_directory": str(folder)})
    if label != "Pregnancy":
        for check in ["inference_verification.json","influence_reconstruction.json"]:
            assert json.loads((folder/check).read_text())["status"] == "passed"
        intervals[label] = pd.read_csv(folder/"adaptive_inference.csv").iloc[0]

summary = pd.DataFrame(decisions)
summary.to_csv(out/"screening_decisions.csv",index=False)
pd.DataFrame(intervals.values()).to_csv(current/"previous_exploratory_inference.csv",index=False)
ra, pregnancy, gnhs = summary.itertuples(index=False)
r, g = intervals["RA"], intervals["GNHS"]
geometry = base / "clinical_quantile_20261105"
checks = json.loads((geometry/"centered_smooth_acceptance.json").read_text())
integration = json.loads((geometry/"centering_acceptance.json").read_text())
text = f"""# 临床分位数估计与患者推断

## 估计目标

CQ_SMOOTH 以临床结局的条件分位数为目标，用中心化核平滑损失估计初始回归。带宽与训练患者数的负五分之一次方成正比。核的第一阶矩为零，在条件密度光滑并具有局部曲率时，总体平滑近似误差为带宽的二阶量。

全部预定临床变量进入可估计设计。高维系数使用 elastic net，随患者数增加而消失的 ridge 项确定数值解。CQ_ORTH 使用原始分位数指示得分，并依据条件密度和工作协方差构造投影。对每位患者汇总 influence 后计算标准误。独立与交换型工作协方差给出的估计按两者在同批患者中的协方差组合。

真实训练资料上的 {checks['checks']} 项独立二次规划检查均通过，目标值最大差异为 {checks['maximum_objective_error']:.2e}，拟合预测的均方根差异最大为 {checks['maximum_prediction_rmse']:.2e}。{integration['cases']} 项确定性总体积分检查核对了核中心化与二阶近似误差。患者交叉拟合的展开、初始估计的乘积条件以及工作协方差的方差界见研究文档。[平滑估计](../../../research/clinical_quantile_20261105/CENTERED_SMOOTHING.md)　[正交展开](../../../research/clinical_quantile_20261105/ORTHOGONAL_EXPANSION.md)　[方差推导](../../../research/clinical_quantile_20261105/THEORY_REVIEW.md)

## 三个队列的预定比较

估计在五个互斥的患者组和三个分位数上进行。各患者的随访损失在患者内部取平均，然后对患者等权平均。预测值按分位数排序并限制为非负。screening_decisions.csv 的 compared_methods 列列出实际完成的方法，authoritative_directory 列给出结果目录。

类风湿关节炎的第 75 百分位用 DAS28-CRP 评分衡量，{r.patients} 名探索患者中，CQ_SMOOTH 的损失为 {ra.new_loss:.6f}，校准 elastic net 为 {ra.best_comparator_loss:.6f}。孕期队列的预定分位数为中位数，以孕周为单位，两个方法的损失分别为 {pregnancy.new_loss:.6f} 和 {pregnancy.best_comparator_loss:.6f}。GNHS 的第 75 百分位用 BMI 衡量，CQ_SMOOTH 的损失为 {gnhs.new_loss:.6f}，临床分位数回归为 {gnhs.best_comparator_loss:.6f}。三个队列均结束于探索比较。[分组损失和比较方法](screening_decisions.csv)

## 患者正交推断

GNHS 的 C_prior_bmi 估计为 {g.estimate:.4f}，患者标准误为 {g.standard_error:.5f}，探索区间为 {g.lower:.4f} 至 {g.upper:.4f}。独立工作协方差给出 {g.independence_standard_error:.5f} 的标准误。使用交换型分量后，估计标准误降低 {g.estimated_se_reduction:.2%}。该系数描述既往 BMI 与随后 BMI 的条件分位数关联。类风湿关节炎队列的 C_age 估计为 {r.estimate:.4f}，患者标准误为 {r.standard_error:.4f}，探索区间为 {r.lower:.4f} 至 {r.upper:.4f}。单名患者占该项估计方差的 {r.largest_patient_variance_share:.2%}，最大五人的合计为 {r.largest_five_patient_variance_share:.2%}。两项研究各自核实了十个投影方向和患者贡献。[逐队列估计与区间](../clinical_followup_20261106/previous_exploratory_inference.csv)

## 结果来源

对孕期队列的 20 个训练设计逐一比较后，入选变量与数值设计一致。本目录的孕期损失仍为该队列的计算来源。类风湿关节炎与 GNHS 的训练设计按患者等权处理缺测及方差，对应的预测和正交推断已经复算。两者的当前比较分别存放在 pw1 目录中。[类风湿关节炎](../cq_ra_westchina_20261105_pw1/discovery_comparisons.csv)　[GNHS](../cq_gnhs_baseline_bmi_20261105_pw1/discovery_comparisons.csv)　[临床随访结果登记](../clinical_followup_20261106/REPORT_zh.md)
"""
(out/"REPORT_zh.md").write_text(text)
