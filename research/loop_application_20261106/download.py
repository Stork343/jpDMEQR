from datetime import datetime, timezone
import json
from pathlib import Path
from zipfile import ZipFile

import boto3
from boto3.s3.transfer import TransferConfig
from botocore import UNSIGNED
from botocore.config import Config

directory = Path("data/raw/loop_20261106")
final = directory / "Loop.zip"
bucket = "live-jchrpublicdatasets"
key = "Diabetes/Public Datasets/Loop study public dataset 2023-01-31.zip"
client = boto3.client("s3", region_name="us-east-1", config=Config(signature_version=UNSIGNED,
    max_pool_connections=36, connect_timeout=20, read_timeout=90, retries={"mode": "standard", "max_attempts": 6},
    proxies={"https": "http://127.0.0.1:7897"}))
size = client.head_object(Bucket=bucket, Key=key)["ContentLength"]
assert size == 1626077053
directory.mkdir(parents=True, exist_ok=True)
if not final.exists():
    incoming = directory / "Loop.s3.incoming"
    client.download_file(bucket, key, str(incoming), Config=TransferConfig(multipart_chunksize=2 * 1024 ** 2,
                         max_concurrency=32, num_download_attempts=8, preferred_transfer_client="classic"))
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
