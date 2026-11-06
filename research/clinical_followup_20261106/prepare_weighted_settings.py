from pathlib import Path
import json

import pandas as pd
import yaml

from panel import registry, settings_directory

audit = Path("results/preflight/clinical_patient_weights_20261106")
new_studies = {v["study_id"]: k for k, v in registry["questions"].items()}
historical = {"iobp2_landmark_20261106": "config/clinical_quantile_iobp2_20261106.yml",
              "loop_night_20261106": "config/clinical_quantile_loop_20261106.yml"}
active = ["ra_westchina_20261105", "pregnancy_detroit_20261105", "gnhs_baseline_bmi_20261105", "gnhs_future_bmi_20261105"]


def write_fixed(path, content):
    text = yaml.safe_dump(content, sort_keys=False)
    if path.exists():
        assert path.read_text() == text
    else:
        path.write_text(text)


for file in sorted(audit.glob("*.json")):
    check = json.loads(file.read_text())
    study_id = check["study_id"]
    if check["changed_panels"] == 0 or study_id.startswith("charite_"):
        continue
    root = Path("results/raw/application") / study_id
    study = yaml.safe_load((root / "implementation/config.yml").read_text())
    if study_id in new_studies:
        method = yaml.safe_load((settings_directory / f"cq_{new_studies[study_id]}.yml").read_text())
    elif study_id in historical:
        method = yaml.safe_load(Path(historical[study_id]).read_text())
    else:
        assert study_id in active
        method = yaml.safe_load(Path("config/clinical_quantile_smooth_20261106.yml").read_text())
    original_report = Path(study["report_directory"])
    report = Path(str(original_report) + "_pw1")
    report.mkdir(parents=True, exist_ok=True)
    if (original_report / "eligibility.json").exists():
        eligibility = json.loads((original_report / "eligibility.json").read_text())
    else:
        rows = pd.read_csv(Path(study["input_directory"]) / "discovery.csv")
        assert rows.role.eq("discovery").all() and rows.groupby("strain_id").outer_fold.nunique().eq(1).all()
        assert set(rows.outer_fold) == set(range(1, study["outer_folds"] + 1))
        eligibility = {"passed": True, "source": "existing frozen discovery dataset and patient folds",
                       "discovery_patients": rows.strain_id.nunique(), "records": len(rows)}
    assert eligibility["passed"]
    if "discovery_patients" not in eligibility or "records" not in eligibility:
        rows = pd.read_csv(Path(study["input_directory"]) / "discovery.csv")
        eligibility.update(discovery_patients=rows.strain_id.nunique(), records=len(rows))
    (report / "eligibility.json").write_text(json.dumps(eligibility, indent=2) + "\n")
    study["output_root"] = str(root) + "_pw1"
    study["report_directory"] = str(report)
    name = "pw_" + study_id
    path = settings_directory / f"{name}.yml"
    write_fixed(path, study)
    method.update(method_id="cq_" + study_id + "_pw1", output_root="results/raw/cq_" + study_id + "_pw1",
                  report_directory="results/preflight/cq_" + study_id + "_pw1", source_studies=[str(path)],
                  inference_quantile=study["decision_quantile"])
    if "inference_target" in study:
        method["inference_target"] = study["inference_target"]
    write_fixed(settings_directory / f"cq_{name}.yml", method)
    print(name)
