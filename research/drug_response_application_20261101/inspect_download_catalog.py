import json
from pathlib import Path

catalog = json.loads(Path("results/scratch/empirical_selection_20261101/cmp_downloads.json").read_text())
records = catalog.get("data", []) + catalog.get("included", [])
for record in records:
    values = record.get("attributes", {})
    if record["type"] != "download_file_group":
        continue
    if values.get("category_id") not in (1, 4, 10):
        continue
    print(json.dumps(values, ensure_ascii=False))
    files = [item["attributes"] for item in records
             if item["type"] == "download_file"
             and item["attributes"]["group_id"] == int(record["id"])]
    for item in sorted(files, key=lambda value: value["timestamp"])[-2:]:
        print(json.dumps(item, ensure_ascii=False))
