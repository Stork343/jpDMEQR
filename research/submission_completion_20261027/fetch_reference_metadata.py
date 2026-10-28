import json
from pathlib import Path
import sys

import requests

doi, key = sys.argv[1:3]
destination = Path("results/scratch/submission_completion_20261027/reference_metadata")
response = requests.get(f"https://api.crossref.org/works/{doi}", timeout=45,
                        headers={"User-Agent": "jpDMEQR manuscript reference audit"})
response.raise_for_status()
document = response.json()
(destination / f"{key}.json").write_text(json.dumps(document, ensure_ascii=False, indent=2))
metadata = document["message"]
print(json.dumps({name: metadata.get(name) for name in
                  ["DOI", "title", "subtitle", "author", "container-title", "volume", "issue", "page", "published-print", "published", "URL"]},
                 ensure_ascii=False, indent=2))
