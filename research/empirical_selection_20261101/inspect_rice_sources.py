import json
from pathlib import Path

import pandas as pd

root = Path("data/raw/rice_baap_20261104")
metadata = json.loads((root / "dataverse_metadata.json").read_text(encoding="utf-8"))
assert metadata["status"] == "OK"
for item in metadata["data"]["latestVersion"]["files"]:
    source = item["dataFile"]
    print(json.dumps({key: source.get(key) for key in ["id", "filename", "filesize", "contentType", "persistentId"]}))
for number in [4, 5]:
    path = root / f"supplementary_{number}.xlsx"
    book = pd.ExcelFile(path)
    for sheet in book.sheet_names:
        data = pd.read_excel(book, sheet_name=sheet, header=None)
        print(path.name, sheet, data.shape)
        print(data.head(7).to_string(index=False, header=False))
