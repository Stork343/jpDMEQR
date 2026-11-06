from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import pandas as pd

base = "https://physionet.org"
release = "/eicu-crd-demo/2.0.1/"
source = Path("data/raw/critical_care_20261106")
source.mkdir(parents=True, exist_ok=True)


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


retrieve("SHA256SUMS.txt")
manifest = pd.read_csv(source / "SHA256SUMS.txt", sep=r"\s+", names=["digest", "filename"])
files = sorted(name for name in manifest.filename if name.endswith(".csv.gz"))
assert len(files) == 31 and {"patient.csv.gz", "lab.csv.gz", "vitalPeriodic.csv.gz"}.issubset(files)
with ThreadPoolExecutor(max_workers=6) as executor:
    sizes = list(executor.map(retrieve, [*files, "LICENSE.txt"]))
sqlite = "sqlite/eicu_v2_0_1.sqlite3.gz"
assert set(manifest.filename) == {*files, "LICENSE.txt", sqlite}
for record in manifest.loc[manifest.filename.ne(sqlite)].itertuples():
    with (source / record.filename).open("rb") as stream:
        assert hashlib.file_digest(stream, "sha256").hexdigest() == record.digest, record.filename
state = {"release": f"{base}/content{release}", "csv_tables": len(files), "compressed_table_bytes": sum(sizes[:-1]),
         "publisher_table_digests_verified": True}
(source / "eicu_release.json").write_text(json.dumps(state, indent=2) + "\n")
print(json.dumps(state), flush=True)
digest, = manifest.loc[manifest.filename.eq(sqlite), "digest"]
subprocess.run([sys.executable, "research/loop_application_20261106/download.py",
                f"{base}/files{release}{sqlite}", str(source / sqlite), "82213572", digest], check=True)
state["publisher_manifest_verified"] = True
(source / "eicu_release.json").write_text(json.dumps(state, indent=2) + "\n")
print(json.dumps(state))
