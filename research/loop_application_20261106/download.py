from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import subprocess
from urllib.parse import quote
from zipfile import ZipFile

directory = Path("data/raw/loop_20261106")
final = directory / "Loop.zip"
bucket = "live-jchrpublicdatasets"
key = "Diabetes/Public Datasets/Loop study public dataset 2023-01-31.zip"
url = f"https://{bucket}.s3.amazonaws.com/{quote(key)}"
size = 1626077053
block = 1024 ** 2
directory.mkdir(parents=True, exist_ok=True)
parts = directory / "ranges"
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


if not final.exists():
    incoming = directory / "Loop.s3.incoming"
    with ThreadPoolExecutor(max_workers=32) as executor:
        chunks = list(executor.map(download_range, range(0, size, block)))
    with incoming.open("wb") as dest:
        for piece in chunks:
            with piece.open("rb") as stream:
                shutil.copyfileobj(stream, dest)
    assert incoming.stat().st_size == size
    with ZipFile(incoming) as archive:
        assert archive.testzip() is None
    incoming.rename(final)
else:
    assert final.stat().st_size == size
    with ZipFile(final) as archive:
        assert archive.testzip() is None
record = {"bucket": bucket, "key": key, "bytes": size, "all_members_verified": True,
          "verified_utc": datetime.now(timezone.utc).isoformat()}
(directory / "verified.json").write_text(json.dumps(record, indent=2) + "\n")
print(json.dumps(record))
