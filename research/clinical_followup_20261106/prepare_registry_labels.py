from pathlib import Path
import yaml

from panel import settings_directory

study = yaml.safe_load((settings_directory/"t1dx_monitoring.yml").read_text())
method = yaml.safe_load((settings_directory/"cq_t1dx_monitoring.yml").read_text())
source = Path(study["input_directory"])
for key in ["input_directory","output_root","report_directory"]:
    study[key] += "_labels1"
inputs = Path(study["input_directory"])
inputs.mkdir(parents=True,exist_ok=True)
report = Path(study["report_directory"])
report.mkdir(parents=True,exist_ok=True)
role = (source/"roles.csv").read_bytes()
if (inputs/"roles.csv").exists():
    assert (inputs/"roles.csv").read_bytes() == role
else:
    (inputs/"roles.csv").write_bytes(role)
settings = settings_directory/"labels_t1dx_monitoring.yml"
method.update(source_studies=[str(settings)])
for key in ["method_id","output_root","report_directory"]:
    method[key] += "_labels1"
for path,data in [(settings,study),(settings_directory/"cq_labels_t1dx_monitoring.yml",method)]:
    text = yaml.safe_dump(data,sort_keys=False)
    if path.exists():
        assert path.read_text() == text
    else:
        path.write_text(text)
