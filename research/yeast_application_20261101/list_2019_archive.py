import json
from pathlib import Path

import requests
from remotezip import RemoteZip

url = "https://www.dropbox.com/sh/jqm7a11zz9laytd/AABaE0EfQxLH6ounPhJ7yYWya?dl=0"
session = requests.Session()
session.proxies.update({"http": "http://127.0.0.1:7897", "https": "http://127.0.0.1:7897"})
with RemoteZip(url, session=session, timeout=90) as archive:
    entries = [{"name": item.filename, "size": item.file_size, "compressed_size": item.compress_size}
               for item in archive.infolist()]
destination = Path("results/preflight/empirical_selection_20261101/yeast_2019_archive.json")
destination.write_text(json.dumps({"source": url, "entries": entries}, indent=2) + "\n", encoding="utf-8")
for entry in entries:
    print(entry["name"], entry["size"])
