import json
import subprocess
import xml.etree.ElementTree as ET
import zipfile
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlparse

from pypdf import PdfReader

root = Path("data/raw/covid_charite_20261105")
metadata = json.loads((root / "PMC8201874.2.json").read_text())
assert metadata["doi"] == "10.1016/j.cels.2021.05.005"
assert metadata["license_code"] == "CC BY" and not metadata["is_manuscript"]
ns = {"s": "http://s3.amazonaws.com/doc/2006-03-01/"}
listing = ET.parse(root / "pmc8201874_listing.xml").getroot()
assert listing.find("s:IsTruncated", ns).text == "false"
sizes = {item.find("s:Key", ns).text: int(item.find("s:Size", ns).text)
         for item in listing.findall("s:Contents", ns)}
sources = [metadata["text_url"], metadata["xml_url"]]
sources += [url for url in metadata["media_urls"] if Path(urlparse(url).path).name.startswith("mmc")]
sources.sort(key=lambda url: Path(urlparse(url).path).name != "mmc7.xlsx")
output = root / "official_cloud"
output.mkdir(exist_ok=True)


def download(source):
    parsed = urlparse(source)
    assert parsed.scheme == "s3" and parsed.netloc == "pmc-oa-opendata"
    key = parsed.path.lstrip("/")
    assert key.startswith("PMC8201874.2/") and key in sizes
    name = Path(key).name
    destination = output / name
    url = "https://pmc-oa-opendata.s3.amazonaws.com/" + key
    if not destination.exists():
        partial = destination.with_name(destination.name + ".partial")
        subprocess.run(["curl", "--silent", "--show-error", "--fail", "--location",
                        "--proxy", "http://127.0.0.1:7897", "--connect-timeout", "20",
                        "--max-time", "180", "--continue-at", "-", "--output", str(partial), url], check=True)
        assert partial.stat().st_size == sizes[key]
        partial.rename(destination)
    assert destination.stat().st_size == sizes[key]
    if destination.suffix in {".xlsx", ".zip"}:
        with zipfile.ZipFile(destination) as archive:
            assert archive.testzip() is None
    if destination.suffix == ".pdf":
        assert len(PdfReader(destination).pages) > 0
    if destination.suffix == ".txt":
        assert metadata["title"] in destination.read_text(encoding="utf-8")
    record = {"file": name, "source": url, "bytes": sizes[key], "status": "verified",
              "source_version": 2, "license": "CC BY", "doi": metadata["doi"],
              "checked_utc": datetime.now(timezone.utc).isoformat()}
    destination.with_name(destination.name + ".source.json").write_text(json.dumps(record, indent=2) + "\n")
    print(f"Verified {name}: {sizes[key]} bytes", flush=True)
    return record


with ThreadPoolExecutor(max_workers=3) as executor:
    records = list(executor.map(download, sources))
manifest = {"status": "downloaded_with_si", "provider": "NIH_NLM_PMC_Cloud_Service",
            "doi": metadata["doi"], "article_version": 2, "files": records}
(output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
print(json.dumps({"status": manifest["status"], "verified_files": len(records)}), flush=True)
