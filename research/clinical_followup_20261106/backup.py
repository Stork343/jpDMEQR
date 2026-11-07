from pathlib import Path
import json
import subprocess

import yaml

from panel import registry, settings_directory

backup = Path("/Users/houjian/jpDMEQR-backups/empirical_selection_20261101/clinical_followup_20261106")
backup.mkdir(parents=True,exist_ok=True)
paths = set()


def completed_files(root):
    for file in root.rglob("*"):
        if not file.is_file() or file.suffix in {".log",".partial",".incoming"} or file.name.startswith("progress"):
            continue
        relative = file.relative_to(root)
        if "models" in relative.parts:
            position = relative.parts.index("models")
            directory = root.joinpath(*relative.parts[:position+2])
            status = directory/"status.json"
            if not status.exists() or json.loads(status.read_text())["status"]!="completed":
                continue
        yield file
for name in ["uwhvf_20261106","critical_care_20261106","gaitrec_20261106","vitaldb_20261106",
             "ispy_20261106","t1dx_20261106","mpc_20261106","commpass_20261106","physionet_sepsis_20261106",
             "vaccine_20261106","trialnet_20261107"]:
    root = Path("data/raw")/name
    for file in root.rglob("*"):
        if not file.is_file() or any(part.endswith(".ranges") or part=="ranges" for part in file.parts):
            continue
        if any(v in file.name for v in [".pypdl", ".incoming",".part",".download.log"]) or file.suffix==".html":
            continue
        if name == "physionet_sepsis_20261106" and file.suffix != ".zip" and file.name not in {"catalog.csv","LICENSE.txt"}:
            continue
        if name == "vaccine_20261106" and file.suffix not in {".rds",".csv",".json",".xml"}:
            continue
        paths.add(file)
for settings in settings_directory.glob("*.yml"):
    paths.add(settings)
    cfg = yaml.safe_load(settings.read_text())
    if "input_directory" not in cfg or "report_directory" not in cfg:
        continue
    report = Path(cfg["report_directory"])
    if not (report/"eligibility.json").exists():
        continue
    paths.update(f for f in Path(cfg["input_directory"]).rglob("*") if f.is_file())
    if not (report/"decision.json").exists():
        continue
    paths.update(completed_files(Path(cfg["output_root"])))
    methodfile = settings_directory/("cq_"+settings.name)
    if not methodfile.exists():
        continue
    method = yaml.safe_load(methodfile.read_text())
    paths.update(completed_files(Path(method["output_root"])))
for location in ["results/raw/clinical_followup_superseded_20261106", "results/preflight/clinical_patient_weights_20261106"]:
    paths.update(f for f in Path(location).rglob("*") if f.is_file())
assert paths
payload = "\n".join(str(f) for f in sorted(paths))+"\n"
subprocess.run(["rsync","-aR","--files-from=-","./",str(backup)+"/"],input=payload,text=True,check=True)
for file in paths:
    target = backup/file
    assert target.is_file() and target.stat().st_size==file.stat().st_size, str(file)
state = {"verified_files":len(paths),"verified_bytes":sum(f.stat().st_size for f in paths),
         "verification":"relative paths and file sizes","clinical_records_published_to_git":False}
out = Path("results/preflight/clinical_followup_20261106")
out.mkdir(parents=True,exist_ok=True)
(out/"backup.json").write_text(json.dumps(state,indent=2)+"\n")
print(json.dumps(state))
