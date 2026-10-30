import argparse
import logging
from pathlib import Path

import aiohttp
from pypdl import Pypdl


parser = argparse.ArgumentParser()
parser.add_argument("url")
parser.add_argument("destination", type=Path)
arguments = parser.parse_args()
arguments.destination.parent.mkdir(parents=True, exist_ok=True)
logging.basicConfig(level=logging.INFO)
downloader = Pypdl(logger=logging.getLogger("GSE65391"))
result = downloader.start(
    arguments.url,
    file_path=str(arguments.destination),
    segments=4,
    retries=0,
    display=False,
    block=True,
    timeout=aiohttp.ClientTimeout(total=7200, sock_connect=30, sock_read=120),
)
if not result or downloader.failed or not arguments.destination.is_file():
    raise RuntimeError("The requested source download did not complete.")
print(f"Downloaded {arguments.destination.name}: {arguments.destination.stat().st_size} bytes")
