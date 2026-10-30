import argparse
import json
import logging
from pathlib import Path
from urllib.request import Request, urlopen

import aiohttp
from pypdl import Pypdl


parser = argparse.ArgumentParser()
parser.add_argument("url")
parser.add_argument("destination", type=Path)
arguments = parser.parse_args()
arguments.destination.parent.mkdir(parents=True, exist_ok=True)
with urlopen(Request(arguments.url, method="HEAD"), timeout=30) as response:
    source_identity = {
        "url": arguments.url,
        "bytes": int(response.headers["Content-Length"]),
        "last_modified": response.headers["Last-Modified"],
    }
if not source_identity["last_modified"]:
    raise RuntimeError("The source response does not identify its modification time.")
identity_path = arguments.destination.with_name(arguments.destination.name + ".source.json")
if identity_path.exists():
    if json.loads(identity_path.read_text()) != source_identity:
        raise RuntimeError("The source file changed during acquisition.")
else:
    identity_path.write_text(json.dumps(source_identity, indent=2))
logging.basicConfig(level=logging.INFO)
downloader = Pypdl(logger=logging.getLogger("GSE65391"))
result = downloader.start(
    arguments.url,
    file_path=str(arguments.destination),
    segments=4,
    retries=0,
    etag_validation=False,
    display=False,
    block=True,
    timeout=aiohttp.ClientTimeout(total=7200, sock_connect=30, sock_read=120),
)
if not result or downloader.failed or not arguments.destination.is_file():
    raise RuntimeError("The requested source download did not complete.")
if arguments.destination.stat().st_size != source_identity["bytes"]:
    raise RuntimeError("The downloaded file size differs from the source response.")
print(f"Downloaded {arguments.destination.name}: {arguments.destination.stat().st_size} bytes")
