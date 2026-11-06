from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from urllib.parse import quote
from zipfile import ZipFile

directory = Path("data/raw/loop_20261106")
final = directory / "Loop.zip"
bucket = "live-jchrpublicdatasets"
key = "Diabetes/Public Datasets/Loop study public dataset 2023-01-31.zip"
url = f"https://{bucket}.s3.amazonaws.com/{quote(key)}"
size = 1626077053
expected_digest = None
if len(sys.argv) > 1:
    if len(sys.argv) not in (4, 5):
        raise ValueError("需要提供来源网址、文件路径、字节数和可选的发布方校验值")
    url, filename, expected_size = sys.argv[1:4]
    final = Path(filename)
    directory = final.parent
    size = int(expected_size)
    if len(sys.argv) == 5:
        expected_digest = sys.argv[4]
        assert expected_digest == "gzip" or len(expected_digest) in (32, 64)
block = 1024 ** 2
directory.mkdir(parents=True, exist_ok=True)
parts = directory / ("ranges" if len(sys.argv) == 1 else final.name + ".ranges")
parts.mkdir(exist_ok=True)


def download_range(start):
    last = min(start + block, size) - 1
    count = last - start + 1
    stored = parts / f"{start:010d}.bin"
    if stored.exists():
        assert stored.stat().st_size == count
        return stored
    incoming = stored.with_suffix(".incoming")
    expected = f"206 bytes {start}-{last}/{size}"
    result = subprocess.run(["curl", "--silent", "--show-error", "--fail", "--location", "--retry", "6",
        "--retry-all-errors", "--retry-delay", "1", "--connect-timeout", "15", "--max-time", "90",
        "--proxy", "http://127.0.0.1:7897", "--range", f"{start}-{last}", "--max-filesize", str(count),
        "--output", str(incoming), "--write-out", "%{response_code} %header{content-range}", url],
        check=True, capture_output=True, text=True)
    assert result.stdout.strip() == expected, (result.stdout, expected)
    assert incoming.stat().st_size == count
    incoming.replace(stored)
    return stored


def verify(path):
    assert path.stat().st_size == size
    if expected_digest is None:
        with ZipFile(path) as archive:
            assert archive.testzip() is None
    elif expected_digest == "gzip":
        import gzip
        with gzip.open(path, "rb") as stream:
            while stream.read(1024 ** 2):
                pass
    else:
        with path.open("rb") as stream:
            algorithm = "md5" if len(expected_digest) == 32 else "sha256"
            assert hashlib.file_digest(stream, algorithm).hexdigest() == expected_digest, str(path)


if not final.exists():
    incoming = final.with_suffix(".s3.incoming")
    with ThreadPoolExecutor(max_workers=int(os.environ.get("CLINICAL_DOWNLOAD_WORKERS", 32))) as executor:
        chunks = list(executor.map(download_range, range(0, size, block)))
    with incoming.open("wb") as dest:
        for piece in chunks:
            with piece.open("rb") as stream:
                shutil.copyfileobj(stream, dest)
    verify(incoming)
    incoming.rename(final)
else:
    verify(final)
record = {"url": url, "bytes": size, "all_members_verified": True,
          "verified_utc": datetime.now(timezone.utc).isoformat()}
if len(sys.argv) == 1:
    record.update(bucket=bucket, key=key)
marker = directory / "verified.json" if len(sys.argv) == 1 else final.with_name(final.name + ".verified.json")
marker.write_text(json.dumps(record, indent=2) + "\n")
print(json.dumps(record))
