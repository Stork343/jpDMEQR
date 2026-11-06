import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import yaml

from panel import registry, settings_directory

root = Path("results/preflight")
out = root / "clinical_followup_20261106"
out.mkdir(exist_ok=True)
question_names = {
    "uwhvf_central": "视野检查", "eicu_lactate": "重症乳酸", "eicu_creatinine": "重症肌酐",
    "gaitrec_propulsion": "患侧步态推进力", "vitaldb_crp": "术后 C 反应蛋白",
    "ispy2_volume": "乳腺癌 MRI 强化体积与基因表达", "ispy2_signaling": "乳腺癌 MRI 强化体积与蛋白信号",
    "t1dx_monitoring": "糖尿病登记随访 HbA1c", "chord_psa": "前列腺癌 PSA", "chord_colorectal": "结直肠癌 CEA",
    "chord_pancreas": "胰腺癌 CA19-9", "chord_breast": "乳腺癌 CA15-3", "commpass_lightchain": "骨髓瘤受累血清轻链",
    "ispy_diameter": "乳腺癌 MRI 最长径", "physionet_ph": "医院 A 血液 pH",
    "vaccine_cxcl10": "流感疫苗接种后的血清 CXCL10"}


def read(path):
    return json.loads(path.read_text())


def verified_comparisons(folder, tau, study_id):
    source = pd.read_csv(folder / "discovery_comparisons.csv")
    x = source.loc[source.study_id.eq(study_id) & source.tau.eq(tau)].set_index("comparator")
    assert not x.index.has_duplicates
    np.testing.assert_allclose(x.relative_improvement, 1-x.current_loss/x.comparator_loss, rtol=0, atol=1e-12)
    assert read(folder / "prediction_acceptance.json")["status"] == "passed"
    return x


records = []
for name, definition in registry["questions"].items():
    study_id = definition["study_id"]
    original = root / "empirical_selection_20261101" / study_id
    roles = read(original / "roles.json")
    assert not roles["reserved_outcomes_used"]
    weighted = settings_directory / f"pw_{study_id}.yml"
    config_name = "labels_t1dx_monitoring" if name == "t1dx_monitoring" else ("pw_"+study_id if weighted.exists() else name)
    study = yaml.safe_load((settings_directory / (config_name+".yml")).read_text())
    method = yaml.safe_load((settings_directory / ("cq_"+config_name+".yml")).read_text())
    report = Path(study["report_directory"])
    row = {"study_id": study_id, "question": question_names[name], "source_doi": definition["source_doi"],
           "decision_quantile": study["decision_quantile"], "reserved_roster": roles["reserved_patients"],
           "report_directory": str(report), "source_roles_directory": str(original),
           "study_settings": str(settings_directory / (config_name+".yml")), "method_directory": method["report_directory"],
           "method_settings": str(settings_directory / ("cq_"+config_name+".yml")),
           "formal_validation_used": False, "manuscript_adoption": False}
    if (report / "eligibility.json").exists():
        q = read(report / "eligibility.json")
        row.update(patients=q["discovery_patients"], records=q["records"], repeated_patients=q["repeated_patients"],
                   phase="awaiting_estimation" if q["passed"] else "closed_at_eligibility")
    else:
        row["phase"] = "source_preparation"
    if (report / "decision.json").exists():
        decision = read(report / "decision.json")
        x = verified_comparisons(Path(method["report_directory"]), study["decision_quantile"], study_id)
        assert read(report / "screen_audit.json")["status"] == "passed"
        assert decision["patients"] == row["patients"] and decision["records"] == row["records"]
        gains = x.loc[["CLINICAL_QR","ENET"], "relative_improvement"]
        np.testing.assert_allclose(gains, [decision["relative_gains"][c] for c in gains.index], atol=1e-12)
        row.update(phase=decision["phase"], gain_over_clinical_qr=gains.CLINICAL_QR, gain_over_enet=gains.ENET,
                   cq_loss=x.current_loss.iloc[0], clinical_qr_loss=x.loc["CLINICAL_QR","comparator_loss"],
                   enet_loss=x.loc["ENET","comparator_loss"])
        if (report / "full_linear_decision.json").exists():
            full = read(report / "full_linear_decision.json")
            np.testing.assert_allclose([full["comparators"][c] for c in x.index], x.relative_improvement, atol=1e-12)
            row.update(phase=full["phase"], full_linear_pass=full["full_linear_pass"],
                       gain_over_joint_direct=x.at["JOINT_DIRECT","relative_improvement"],
                       gain_over_joint_offset=x.at["JOINT_OFFSET","relative_improvement"],
                       best_checked_comparator=x.comparator_loss.idxmin())
        if (report / "clinical_enet.json").exists():
            extra = read(report / "clinical_enet.json")
            np.testing.assert_allclose(extra["gain"], 1-row["cq_loss"]/extra["loss"], atol=1e-12)
            row["gain_over_clinical_enet"] = extra["gain"]
    records.append(row)
frame = pd.DataFrame(records)
frame.to_csv(out / "screening_register.csv", index=False)
fitted = frame.dropna(subset=["cq_loss"]).copy()
largest = fitted.assign(minimum=fitted[["gain_over_clinical_qr","gain_over_enet"]].min(axis=1)).sort_values("minimum").iloc[-1]
pass_count = int((fitted[["gain_over_clinical_qr","gain_over_enet"]].min(axis=1) >= registry["required_relative_gain"]).sum())
document = [
    "# 临床随访比较结果", "",
    f"本次登记 {len(frame)} 个问题，{len(fitted)} 个已经完成预定分位数的模型比较，"
    f"{frame.phase.eq('closed_at_eligibility').sum()} 个结束于资料质量检查。"
    f"在两项必要比较中均达到 2% 改善的问题有 {pass_count} 个。患者数按问题计算；部分问题来自同一批患者，不能把人数相加。", "",
    "探索使用按患者分组的五次交叉验证。每名患者在各次随访的损失先取平均，再对患者取平均。正数表示 CQ_SMOOTH 的损失较低。"
    "临床分位数回归使用预定临床资料和目标变量，校准 elastic net 还使用相应的高维资料。"
    "两个必要比较只是完整比较的第一步。", "",
    f"两项增幅中较小者最高的是{largest.question}，相对临床分位数回归改善 {largest.gain_over_clinical_qr:.2%}，"
    f"相对校准 elastic net 改善 {largest.gain_over_enet:.2%}。各题的患者数、分位数、两种增幅、绝对损失和来源 DOI 见"
    " [screening_register.csv](screening_register.csv)。", "",
    "## 数值复核", "",
    "按患者计算变量观测率与训练方差后，逐一复核了 224 个历史训练设计。发生变化且有探索比较的六个问题均已重新拟合。"
    "I-SPY 2 基因表达与 CoMMpass 均在 20 个训练设计中独立重选了登记的 300 个基因。"
    "糖尿病登记资料按发布方 F、M、T 编码构建性别指标并重新拟合。独立计算核实了 HbA1c 和监测次数差。", "",
    "结直肠癌、胰腺癌、乳腺癌及前列腺癌的标志物均通过独立程序从原始化验与肿瘤测序记录重建，"
    "同时逐人核对了取样日期和基因交互项。骨髓瘤的轻链由 GDC 的逐次检查独立重建，另从 Xena 表达矩阵重新提取五个增殖基因。", "",
    "## 研究状态", "",
    "登记保留的随访数值仍封存。验证尝试表中已有的 Charité 记录属于此前的研究；本轮没有使用它的验证患者。"
    "临床发现与方法优越性须由预定比较以及独立患者区间支持。[验证登记](../empirical_selection_20261101/validation_attempts.csv)", "",
    "患者划分与数值结果对应的目录在结果表的 report_directory、source_roles_directory 两列。"
    "带有 pw1 的目录采用核实后的患者权重，labels1 记录糖尿病登记资料的复算。完整患者记录和原始拟合保存在本机。"]
vaccine = frame.loc[frame.study_id.eq("vaccine_cxcl10_20261106")]
if len(vaccine) and pd.notna(vaccine.iloc[0].get("gain_over_joint_direct")):
    trial = vaccine.iloc[0]
    inference = pd.read_csv(Path(trial.method_directory)/"adaptive_inference.csv").iloc[0]
    assert read(Path(trial.method_directory)/"inference_verification.json")["status"] == "passed"
    assert read(Path(trial.method_directory)/"influence_reconstruction.json")["status"] == "passed"
    document.extend(["", "## 流感疫苗接种后的血清 CXCL10", "",
        f"两项接种研究共提供 {int(trial.patients)} 名探索参与者的 {int(trial.records)} 次配对测量。"
        f"CQ_SMOOTH 相对临床分位数回归改善 {trial.gain_over_clinical_qr:.2%}，相对全变量校准 elastic net 改善 "
        f"{trial.gain_over_enet:.2%}。JOINT_DIRECT 的损失低于 CQ_SMOOTH；按 JOINT_DIRECT 的损失作分母，"
        f"CQ_SMOOTH 高出 {-trial.gain_over_joint_direct:.2%}。JOINT_OFFSET 对应的差值是 "
        f"{-trial.gain_over_joint_offset:.2%}。完整线性比较据此结束。"
        " [线性比较](../empirical_selection_20261101/vaccine_cxcl10_20261106/full_linear_decision.json)", "",
        "预定临床系数描述同一参与者干扰素表达指标较接种前升高十个百分点时，与当次血清浓度条件分位数的关联。"
        f"探索估计为 {inference.estimate:.3f}，患者标准误为 {inference.standard_error:.3f}，"
        f"区间为 {inference.lower:.3f} 至 {inference.upper:.3f}，跨过零。"
        f"{int(trial.reserved_roster)} 名验证参与者保持封存。"
        " [患者推断](../cq_vaccine_cxcl10_20261106/adaptive_inference.csv)"])
(out / "REPORT_zh.md").write_text("\n".join(document)+"\n")

subprocess.run([sys.executable, "research/clinical_quantile_20261105/write_smooth_report.py"], check=True)
