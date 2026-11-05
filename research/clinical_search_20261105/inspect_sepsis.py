import json
import zipfile
from pathlib import Path

import pandas as pd

root = Path("data/raw/sepsis_nrw_20261105")
records = []
for name in ["basic_patient_information.xlsx", "full_protein_information.xlsx", "normalized_protein_intensities.xlsx"]:
    destination = root / name
    partial = destination.with_name(destination.name + ".partial")
    if not destination.exists():
        with zipfile.ZipFile(partial) as archive:
            assert archive.testzip() is None
        partial.rename(destination)
    workbook = pd.ExcelFile(destination)
    records.append({"file": name, "bytes": destination.stat().st_size, "status": "verified",
                    "source": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2024/12/PXD055932/" + name})
    print(json.dumps({"file": name, "sheets": workbook.sheet_names}))
    for sheet in workbook.sheet_names:
        frame = pd.read_excel(workbook, sheet_name=sheet)
        print(json.dumps({"sheet": sheet, "shape": frame.shape, "columns": [str(x) for x in frame.columns[:25]]}))
        print(frame.iloc[:5, :8].to_string(index=False))
(root / "pride_manifest.json").write_text(json.dumps(records, indent=2) + "\n")
