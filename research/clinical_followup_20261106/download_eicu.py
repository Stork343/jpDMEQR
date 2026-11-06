from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess

from bs4 import BeautifulSoup
import requests

base = "https://physionet.org"
release = "/eicu-crd-demo/2.0.1/"
source = Path("data/raw/critical_care_20261106")
source.mkdir(parents=True, exist_ok=True)
response = requests.get(f"{base}/content{release}", timeout=40, proxies={"https": "http://127.0.0.1:7897"})
response.raise_for_status()
page = BeautifulSoup(response.content, "html.parser")
files = sorted({a.get_text(strip=True) for a in page.select("a[href]") if a.get_text(strip=True).endswith(".csv.gz")})
assert len(files) >= 30 and {"patient.csv.gz", "lab.csv.gz", "vitalPeriodic.csv.gz"}.issubset(files)


def retrieve(filename):
    path = source / filename
    if not path.exists():
        partial = path.with_name(path.name + ".part")
        subprocess.run(["curl", "--silent", "--show-error", "--fail", "--location", "--retry", "5",
                        "--retry-all-errors", "--connect-timeout", "18", "--retry-delay", "1", "--max-time", "180",
                        "--proxy", "http://127.0.0.1:7897", "--continue-at", "-", "--output", str(partial),
                        f"{base}/files{release}{filename}"], check=True)
        partial.rename(path)
    return path.stat().st_size


with ThreadPoolExecutor(max_workers=6) as executor:
    sizes = list(executor.map(retrieve, [*files, "LICENSE.txt", "SHA256SUMS.txt"]))
subprocess.run(["shasum", "-a", "256", "--check", "--status", "SHA256SUMS.txt"], cwd=source, check=True)
state = {"release": f"{base}/content{release}", "csv_tables": len(files), "compressed_table_bytes": sum(sizes[:-2]),
         "publisher_manifest_verified": True}
(source / "eicu_release.json").write_text(json.dumps(state, indent=2) + "\n")
print(json.dumps(state))
