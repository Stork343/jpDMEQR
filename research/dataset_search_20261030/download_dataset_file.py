import argparse
import json
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from zipfile import ZipFile


parser = argparse.ArgumentParser()
parser.add_argument("url")
parser.add_argument("destination", type=Path)
parser.add_argument("expected_bytes", type=int)
arguments = parser.parse_args()
destination = arguments.destination.resolve()
if destination.exists():
    raise FileExistsError(destination)
destination.parent.mkdir(parents=True, exist_ok=True)
partial = destination.with_name(destination.name + ".downloading")
parts = destination.parent / (destination.name + ".parts")
parts.mkdir(exist_ok=True)
identity = {"url": arguments.url, "bytes": arguments.expected_bytes}
manifest = parts / "source.json"
if manifest.exists() and json.loads(manifest.read_text()) != identity:
    raise ValueError("The saved transfer belongs to a different source.")
manifest.write_text(json.dumps(identity, indent=2))
chunk_size = 1024 * 1024
ranges = [(start, min(start + chunk_size, arguments.expected_bytes) - 1)
          for start in range(0, arguments.expected_bytes, chunk_size)]


def transfer(item):
    index, (start, end) = item
    path = parts / f"part_{index:04d}"
    if path.exists() and path.stat().st_size == end - start + 1:
        return path
    subprocess.run([
        "curl", "--fail", "--location", "--silent", "--show-error",
        "--connect-timeout", "20", "--max-time", "180", "--retry", "3",
        "--retry-all-errors", "--range", f"{start}-{end}",
        "--output", str(path), arguments.url,
    ], check=True)
    if path.stat().st_size != end - start + 1:
        raise ValueError(f"Incomplete source transfer: {index}")
    print(f"Verified segment {index + 1}/{len(ranges)}", flush=True)
    return path


with ThreadPoolExecutor(max_workers=4) as executor:
    completed = list(executor.map(transfer, enumerate(ranges)))
with partial.open("wb") as output:
    for path in completed:
        with path.open("rb") as source:
            shutil.copyfileobj(source, output)
if partial.stat().st_size != arguments.expected_bytes:
    raise ValueError("The file size differs from the repository metadata.")
if destination.suffix.lower() in (".zip", ".xlsx"):
    with ZipFile(partial) as archive:
        if archive.testzip() is not None:
            raise ValueError("The downloaded archive is incomplete.")
partial.rename(destination)
print(f"Verified {destination.name}: {destination.stat().st_size} bytes")
