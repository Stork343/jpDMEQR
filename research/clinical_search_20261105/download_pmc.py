import argparse
import json
import subprocess
import xml.etree.ElementTree as ET
import zipfile
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

from pypdf import PdfReader

parser = argparse.ArgumentParser()
parser.add_argument("--root", type=Path, required=True)
parser.add_argument("--article", required=True)
parser.add_argument("--doi", required=True)
parser.add_argument("--si", action="store_true", required=True)
args = parser.parse_args()
root = args.root
listing = ET.parse(root / "pmc_listing.xml").getroot()
namespace = {"s": "http://s3.amazonaws.com/doc/2006-03-01/"}
assert listing.find("s:IsTruncated", namespace).text == "false"
objects = {entry.find("s:Key", namespace).text: int(entry.find("s:Size", namespace).text)
           for entry in listing.findall("s:Contents", namespace)}
prefix = args.article + "/"
objects = {key: size for key, size in objects.items() if key.startswith(prefix)}
assert objects
output = root / "official_cloud"
output.mkdir(exist_ok=True)


def fetch(key):
    destination = output / Path(key).name
    url = "https://pmc-oa-opendata.s3.us-east-1.amazonaws.com/" + key
    if not destination.exists():
        partial = destination.with_name(destination.name + ".partial")
        if partial.exists():
            assert partial.stat().st_size <= objects[key]
        if not partial.exists() or partial.stat().st_size < objects[key]:
            subprocess.run(["curl", "--silent", "--show-error", "--fail", "--location",
                            "--proxy", "http://127.0.0.1:7897", "--connect-timeout", "20",
                            "--max-time", "600", "--continue-at", "-", "--output", str(partial), url], check=True)
        assert partial.stat().st_size == objects[key]
        partial.rename(destination)
    assert destination.stat().st_size == objects[key]
    if destination.suffix in {".xlsx", ".zip", ".docx"}:
        with zipfile.ZipFile(destination) as archive:
            assert archive.testzip() is None
    if destination.suffix == ".pdf":
        reader = PdfReader(destination)
        assert len(reader.pages) > 0
    if destination.suffix == ".xml":
        ET.parse(destination)
    record = {"file": destination.name, "source": url, "bytes": objects[key], "status": "verified",
              "article": args.article, "doi": args.doi, "checked_utc": datetime.now(timezone.utc).isoformat()}
    destination.with_name(destination.name + ".source.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps({"file": destination.name, "status": "verified", "bytes": objects[key]}), flush=True)
    return record


metadata_key = prefix + args.article + ".json"
metadata_record = fetch(metadata_key)
metadata = json.loads((output / Path(metadata_key).name).read_text())
assert metadata["doi"] == args.doi and not metadata["is_manuscript"]
assert not metadata["is_retracted"]
keys = [key for key in objects if Path(key).suffix in {".txt", ".xml", ".pdf", ".xlsx", ".xls", ".csv", ".tsv", ".zip", ".docx"}]
keys.sort(key=lambda key: (Path(key).suffix not in {".txt", ".xlsx"}, objects[key]))
with ThreadPoolExecutor(max_workers=3) as executor:
    records = list(executor.map(fetch, keys))
manifest = {"status": "downloaded_with_si", "article": args.article, "doi": args.doi,
            "license": metadata["license_code"], "files": [metadata_record] + records}
(output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
print(json.dumps({"status": manifest["status"], "files": len(manifest["files"])}), flush=True)
