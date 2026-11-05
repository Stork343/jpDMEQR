import argparse
import json
import subprocess
import zipfile
from pathlib import Path
from urllib.parse import urljoin, urlparse

import pandas as pd
from bs4 import BeautifulSoup

parser = argparse.ArgumentParser()
parser.add_argument("tables", nargs="+", choices=[f"tableS{i}" for i in range(1, 8)])
args = parser.parse_args()
root = Path("data/raw/gnhs_20261105")
page = BeautifulSoup((root / "resource_page.html").read_text(), "html.parser")
links = {link.get_text(strip=True): urljoin("https://guomics.com/publications/2023crm_caix/", link["href"])
         for link in page.find_all("a", href=True) if link.get_text(strip=True) in args.tables}
assert set(links) == set(args.tables)
for label in args.tables:
    url = links[label]
    assert url.startswith("https://guomics-web.oss-cn-beijing.aliyuncs.com/")
    extension = Path(urlparse(url).path).suffix
    assert extension == ".zip"
    destination = root / f"{label}{extension}"
    if not destination.exists():
        partial = destination.with_name(destination.name + ".partial")
        subprocess.run(["curl", "--silent", "--show-error", "--fail", "--location",
                        "--proxy", "http://127.0.0.1:7897", "--connect-timeout", "20",
                        "--max-time", "600", "--continue-at", "-", "--output", str(partial), url], check=True)
        with zipfile.ZipFile(partial) as archive:
            assert archive.testzip() is None
        partial.rename(destination)
    with zipfile.ZipFile(destination) as archive:
        assert archive.testzip() is None
    metadata = {"file": destination.name, "source": url, "bytes": destination.stat().st_size,
                "status": "verified", "source_doi": "10.1016/j.xcrm.2023.101172"}
    (root / f"{label}.source.json").write_text(json.dumps(metadata, indent=2) + "\n")
    with zipfile.ZipFile(destination) as archive:
        members = [item for item in archive.infolist() if not item.is_dir()]
        assert members and all(Path(item.filename).suffix == ".xlsx" for item in members)
        for item in members:
            assert not Path(item.filename).is_absolute() and ".." not in Path(item.filename).parts
            directory = root / "tables" / label
            directory.mkdir(parents=True, exist_ok=True)
            path = Path(archive.extract(item, directory))
            workbook = pd.ExcelFile(path)
            print(json.dumps({"file": path.name, "sheets": workbook.sheet_names}), flush=True)
            for sheet in workbook.sheet_names:
                frame = pd.read_excel(workbook, sheet_name=sheet)
                print(json.dumps({"sheet": sheet, "shape": frame.shape, "columns": [str(x) for x in frame.columns[:25]]}), flush=True)
                print(frame.iloc[:4, :8].to_string(index=False), flush=True)
