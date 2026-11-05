import json
import zipfile
from pathlib import Path

import pandas as pd

root = Path("data/raw/ra_westchina_20261105/official_cloud")
for path in sorted(root.glob("*.xlsx")):
    workbook = pd.ExcelFile(path)
    print(json.dumps({"file": path.name, "sheets": workbook.sheet_names}))
    for sheet in workbook.sheet_names:
        frame = pd.read_excel(workbook, sheet_name=sheet)
        print(json.dumps({"sheet": sheet, "shape": frame.shape,
                          "columns": [str(x) for x in frame.columns[:18]],
                          "last_columns": [str(x) for x in frame.columns[-8:]]}))
        print(frame.iloc[:4, :min(8, len(frame.columns))].to_string(index=False))
with zipfile.ZipFile(root / "41467_2025_62032_MOESM7_ESM.zip") as archive:
    print(json.dumps({"source_archive": [item.filename for item in archive.infolist() if not item.is_dir()]}))
clinical = pd.read_excel(root / "41467_2025_62032_MOESM3_ESM.xlsx", header=2)
print(clinical.head(3).to_string(index=False))
for column in clinical.columns:
    values = clinical[column]
    print(json.dumps({"column": column, "missing": int(values.isna().sum()),
                      "unique": int(values.nunique()), "examples": values.dropna().astype(str).unique()[:12].tolist()}))
print(clinical.groupby("Group").agg(records=("Sample ID", "size"), patients=("Patient ID", "nunique")).to_string())
print(clinical.loc[clinical["Patient ID"].duplicated(keep=False)].head(15).to_string(index=False))
