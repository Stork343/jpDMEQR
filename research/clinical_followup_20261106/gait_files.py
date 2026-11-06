from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import subprocess
import sys

from panel import study_config

study = study_config("gaitrec_propulsion")
root = Path(study["source_directory"])
root.mkdir(parents=True, exist_ok=True)
catalog = root / "index"
catalog.mkdir(exist_ok=True)


def public_json(url, destination):
    if not destination.exists():
        partial = destination.with_suffix(".incoming")
        subprocess.run(["curl", "--silent", "--show-error", "--fail", "--location", "--retry", "6",
                        "--retry-all-errors", "--retry-delay", "1", "--connect-timeout", "15", "--max-time", "40",
                        "--proxy", "http://127.0.0.1:7897", "--output", str(partial), url], check=True)
        json.loads(partial.read_text())
        partial.rename(destination)
    return json.loads(destination.read_text())


collection = public_json("https://api.figshare.com/v2/collections/4788012/articles?page_size=100",
                         root / "collection.json")


def file_record(article):
    details = public_json(f"https://api.figshare.com/v2/articles/{article['id']}", catalog / f"{article['id']}.json")
    assert details["license"]["name"] == "CC0"
    file, = details["files"]
    return {"filename": file["name"], "size": file["size"], "url": file["download_url"],
            "md5": file["supplied_md5"], "article": article["id"], "doi": details["doi"]}


with ThreadPoolExecutor(max_workers=5) as executor:
    files = list(executor.map(file_record, [a for a in collection if a["title"].startswith("GRF_")]))
assert len(files) == 21
(root / "files.json").write_text(json.dumps(files, indent=2) + "\n")
mode = sys.argv[1]
if mode == "index":
    print(json.dumps({key: sum(f["size"] for f in files if key in f["filename"]) for key in ["_PRO_", "_RAW_", "metadata"]}))
    raise SystemExit(0)
assert mode in {"processed", "raw"}
selected = [f for f in files if ("_PRO_" if mode == "processed" else "_RAW_") in f["filename"]]
assert len(selected) == 10
environment = dict(os.environ)
environment.setdefault("CLINICAL_DOWNLOAD_WORKERS", "10")


def retrieve(file):
    path = root / file["filename"]
    subprocess.run([sys.executable, "research/loop_application_20261106/download.py", file["url"], str(path),
                    str(file["size"]), file["md5"]], env=environment, check=True)
    print(file["filename"], flush=True)


with ThreadPoolExecutor(max_workers=int(os.environ.get("GAIT_FILE_WORKERS", 3))) as executor:
    list(executor.map(retrieve, selected))
