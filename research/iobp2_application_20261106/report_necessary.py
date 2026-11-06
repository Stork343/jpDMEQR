import json
from pathlib import Path
import sys

import pandas as pd
import yaml

method = yaml.safe_load(Path(sys.argv[1]).read_text())
study = yaml.safe_load(Path(method["source_studies"][0]).read_text())
report = Path(study["report_directory"])
output = Path(method["report_directory"])
roster = json.loads((report / "roles.json").read_text())
panel = json.loads((report / "panel_verification.json").read_text())
base_check = json.loads((report / "screen_audit.json").read_text())
model_check = json.loads((output / "prediction_acceptance.json").read_text())
assert base_check["status"] == model_check["status"] == "passed" and base_check["models_checked"] == 10
comparisons = pd.read_csv(output / "discovery_comparisons.csv").query("decision_quantile").set_index("comparator")
assert set(comparisons.index) == {"CLINICAL_QR", "ENET"}
gain = comparisons.relative_improvement
passed = bool(gain.min() >= study["screening_minimum_relative_gain"])
diagnostics = pd.read_csv(output / f"{study['study_id']}_fit_diagnostics.csv")
chosen = diagnostics.loc[diagnostics.tau.eq(study["decision_quantile"]), ["fold", "alpha", "selected_index",
                         "nonzero_proteins", "clinical_columns", "maximum_kkt"]]
chosen = chosen.rename(columns={"nonzero_proteins": "nonzero_glucose_profile_features"})
chosen.to_csv(output / "decision_folds.csv", index=False)
counts = chosen.nonzero_glucose_profile_features.to_list()
state = {"study_id": study["study_id"], "discovery_patients": panel["discovery_patients_reconstructed"],
         "discovery_measurements": panel["records_reconstructed"], "reserved_roster": roster["reserved_patients"],
         "decision_quantile": study["decision_quantile"], "relative_gains": gain.to_dict(),
         "profile_features_per_fold": counts, "necessary_comparisons_pass": passed,
         "phase": "eligible_for_additional_comparisons" if passed else "closed_after_necessary_comparisons",
         "formal_validation_used": False, "manuscript_adoption": False}
(report / "decision.json").write_text(json.dumps(state, indent=2) + "\n")
is_loop = study.get("single_arm_study", False)
name = "Loop" if is_loop else "IOBP2"
attribution = ("The source of the data is the Loop Study (sponsored by the Jaeb Center for Health Research and "
               "funded by the Helmsley Charitable Trust), but the analyses, content and conclusions presented "
               "herein are solely the responsibility of the authors and have not been reviewed or approved by the study sponsor."
               if is_loop else "The source of the data is the Insulin Only Bionic Pancreas Pivotal Trial, but the analyses, "
               "content and conclusions presented herein are solely the responsibility of the authors and have not been "
               "reviewed or approved by the Bionic Pancreas Research Group or Beta Bionics.")
source = ("Lum 等，Diabetes Technology & Therapeutics，2021，DOI 10.1089/dia.2020.0535。" if is_loop
          else "Bionic Pancreas Research Group，New England Journal of Medicine，2022，DOI 10.1056/NEJMoa2205225。")
lines = [f"# {name} 入组两周时的夜间血糖预测", "",
         f"探索侧登记 {roster['discovery_patients']} 人，达到监测条件的患者有 {state['discovery_patients']} 人，"
         f"每人提供三个时段，共 {state['discovery_measurements']} 个结局。{roster['reserved_patients']} 名患者留作验证。",
         "", "到第 14 日结束时，以已经取得的监测数据预测随后三个 21 日时段的夜间高值。"
         "每夜计算最高 10% 血糖读数的均值，每个时段对合格夜间等权平均。评价以患者等权的第 75 百分位预测损失为准，数值越小越好。",
         "", f"CQ_SMOOTH 的损失为 {comparisons.current_loss.iloc[0]:.4f}。"
         f"临床分位数回归为 {comparisons.loc['CLINICAL_QR','comparator_loss']:.4f}，"
         f"经分位数校准的 elastic net 为 {comparisons.loc['ENET','comparator_loss']:.4f}。"
         f"相对改善分别为 {gain['CLINICAL_QR']:+.2%} 和 {gain['ENET']:+.2%}。",
         "", f"每折训练结束后，288 个血糖曲线变量中系数绝对值超过 10⁻⁷ 的个数依次为 {', '.join(map(str,counts))}。"
         "Python 与 R 分别重建了预测误差，DuckDB 从五分钟时格重新计算了夜间结局和 288 个输入变量。",
         "", ("两项必要比较均达到 2% 条件，尚须完成登记的其余比较。" if passed else
              "必要比较未达到登记的 2% 条件。保留患者的随访结局继续封存，本次计算未形成经独立验证的临床发现。"),
         "", source, "", attribution, ""]
(output / "REPORT_zh.md").write_text("\n".join(lines))
print(json.dumps(state))
