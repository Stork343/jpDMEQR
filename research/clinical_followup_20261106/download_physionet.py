from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from threading import local
from zipfile import ZIP_DEFLATED, ZipFile
import sys

from lxml import etree
import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

root = Path("data/raw/physionet_sepsis_20261106")
root.mkdir(parents=True, exist_ok=True)
base = "https://physionet-open.s3.amazonaws.com"
cache = local()


def session():
    if not hasattr(cache, "client"):
        cache.client = requests.Session()
        retry = Retry(total=4, backoff_factor=.8, status_forcelist=[429, 500, 502, 503, 504])
        cache.client.mount("https://", HTTPAdapter(max_retries=retry))
    return cache.client


def inventory(label):
    prefix = "challenge-2019/1.0.0/training/training_set" + label + "/"
    params = {"list-type": "2", "prefix": prefix}
    records = []
    while True:
        response = session().get(base, params=params, timeout=(12, 50))
        response.raise_for_status()
        page = etree.fromstring(response.content)
        for item in page.findall("{*}Contents"):
            key = item.findtext("{*}Key")
            if key.endswith(".psv"):
                records.append((key, int(item.findtext("{*}Size")), label, Path(key).stem))
        token = page.findtext("{*}NextContinuationToken")
        if token is None:
            break
        params["continuation-token"] = token
    return pd.DataFrame(records, columns=["key", "bytes", "site", "patient"]).sort_values("patient")


def retrieve(row):
    response = session().get(base + "/" + row.key, timeout=(12, 45))
    response.raise_for_status()
    assert len(response.content) == row.bytes
    return row.patient + ".psv", response.content


catalog_file = root / "catalog.csv"
if catalog_file.exists():
    catalog = pd.read_csv(catalog_file)
else:
    with ThreadPoolExecutor(max_workers=2) as pool:
        catalog = pd.concat(pool.map(inventory, ["A", "B"]), ignore_index=True)
    assert catalog.patient.is_unique and catalog.groupby("site").size().to_dict() == {"A": 20336, "B": 20000}
    catalog["batch"] = catalog.groupby("site").cumcount() // 250
    incoming = root / "catalog.incoming"
    catalog.to_csv(incoming, index=False)
    incoming.replace(catalog_file)
print("Public roster", catalog.groupby("site").size().to_dict(), flush=True)
if sys.argv[1:] == ["inventory"]:
    raise SystemExit(0)
requested = sys.argv[1:] or ["A", "B"]
assert set(requested) <= {"A", "B"}
with ThreadPoolExecutor(max_workers=24) as pool:
    for (site, batch), group in catalog.loc[catalog.site.isin(requested)].groupby(["site", "batch"]):
        final = root / f"{site}_{batch:03}.zip"
        if not final.exists():
            incoming = final.with_suffix(".incoming")
            with ZipFile(incoming, "w", compression=ZIP_DEFLATED, compresslevel=3) as archive:
                for name, content in pool.map(retrieve, group.itertuples(index=False)):
                    archive.writestr(name, content)
            incoming.replace(final)
        with ZipFile(final) as archive:
            sizes = {entry.filename: entry.file_size for entry in archive.infolist()}
        assert sizes == dict(zip(group.patient + ".psv", group.bytes))
        print(site, batch, len(group), flush=True)
