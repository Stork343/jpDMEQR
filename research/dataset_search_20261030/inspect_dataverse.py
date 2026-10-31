import argparse
import json
from pathlib import Path


parser = argparse.ArgumentParser()
parser.add_argument("metadata", type=Path)
arguments = parser.parse_args()
record = json.loads(arguments.metadata.read_text())
version = record["data"]["latestVersion"]
fields = version["metadataBlocks"]["citation"]["fields"]
print(json.dumps({field["typeName"]: field["value"] for field in fields}, ensure_ascii=False, indent=2))
print(json.dumps({key: version.get(key) for key in ("license", "termsOfUse", "termsOfAccess", "restrictions")}, ensure_ascii=False, indent=2))
print(json.dumps([
    {"id": item["dataFile"]["id"], "name": item["dataFile"]["filename"],
     "size": item["dataFile"]["filesize"], "restricted": item.get("restricted"),
     "description": item["dataFile"].get("description")}
    for item in version["files"]
], ensure_ascii=False, indent=2))
