import argparse
import json
import logging
import os
from datetime import datetime, timezone
from pathlib import Path

from aiohttp import ClientTimeout

parser = argparse.ArgumentParser()
parser.add_argument("url")
parser.add_argument("destination", type=Path)
parser.add_argument("--bytes", type=int)
parser.add_argument("--segments", type=int, default=8)
parser.add_argument("--proxy", default="http://127.0.0.1:7897")
args = parser.parse_args()
destination = args.destination.resolve()
if destination.exists():
    raise FileExistsError(destination)
destination.parent.mkdir(parents=True, exist_ok=True)
runtime_dir = Path(__file__).resolve().parents[2] / "results/scratch/runtime"
runtime_dir.mkdir(parents=True, exist_ok=True)
os.chdir(runtime_dir)
from pypdl import Pypdl

partial = destination.with_name(destination.name + ".pypdl.part")
logger = logging.getLogger(str(destination))
logger.setLevel(logging.WARNING)
logger.addHandler(logging.FileHandler(destination.with_name(destination.name + ".download.log")))
downloader = Pypdl(logger=logger)
result = downloader.start(
    args.url,
    file_path=str(partial),
    segments=args.segments,
    retries=0,
    display=False,
    proxy=args.proxy,
    timeout=ClientTimeout(total=None, connect=45, sock_read=180),
)
if not result or not downloader.completed or downloader.failed:
    raise RuntimeError("数据下载未完成")
size = partial.stat().st_size
if args.bytes is not None and size != args.bytes:
    raise ValueError(f"文件长度不符合要求: {size} != {args.bytes}")
partial.rename(destination)
record = {
    "url": args.url,
    "path": str(destination),
    "bytes": size,
    "downloaded_utc": datetime.now(timezone.utc).isoformat(),
    "completed": True,
}
destination.with_name(destination.name + ".download.json").write_text(
    json.dumps(record, ensure_ascii=False, indent=2), encoding="utf-8"
)
print(json.dumps(record, ensure_ascii=False))
