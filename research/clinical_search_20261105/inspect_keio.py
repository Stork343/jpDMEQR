import json
from pathlib import Path

import pandas as pd

root = Path("data/raw/ra_keio_20261105/official_cloud")
for path in sorted(root.glob("*.xlsx")):
    workbook = pd.ExcelFile(path)
    print(json.dumps({"file": path.name, "sheets": workbook.sheet_names}))
    for sheet in workbook.sheet_names:
        frame = pd.read_excel(workbook, sheet_name=sheet)
        print(json.dumps({"sheet": sheet, "shape": frame.shape, "columns": [str(x) for x in frame.columns[:12]]}))
        print(frame.iloc[:3, :6].to_string(index=False))
