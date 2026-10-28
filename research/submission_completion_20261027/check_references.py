import csv
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import unicodedata

import bibtexparser
import requests
from pylatexenc.latex2text import LatexNodes2Text
from lxml import html

root = Path.cwd()
destination = root / "results/scratch/submission_completion_20261027/reference_metadata"
destination.mkdir(parents=True, exist_ok=True)
bibliography = bibtexparser.parse_file(root / "manuscript/references.bib")
if bibliography.failed_blocks:
    raise ValueError("The reference file contains an unparsed BibTeX block")


def field(entry, name):
    item = entry.fields_dict.get(name)
    return item.value if item is not None else ""


def normalize(value):
    plain = html.fragment_fromstring(value, create_parent=True).text_content() if value else ""
    text = LatexNodes2Text().latex_to_text(plain)
    text = unicodedata.normalize("NFKD", text).casefold()
    return "".join(character for character in text if character.isalnum())


def fetch(entry):
    doi = LatexNodes2Text().latex_to_text(field(entry, "doi"))
    path = destination / f"{entry.key}.json"
    if path.exists():
        document = json.loads(path.read_text())
    else:
        response = requests.get(f"https://api.crossref.org/works/{doi}", timeout=45,
                                headers={"User-Agent": "jpDMEQR manuscript reference audit"})
        response.raise_for_status()
        document = response.json()
        path.write_text(json.dumps(document, ensure_ascii=False, indent=2))
    metadata = document["message"]
    if metadata["DOI"].casefold() != doi.casefold():
        raise ValueError(f"The returned DOI does not match {entry.key}")
    date = metadata.get("published-print", metadata.get("published", metadata.get("issued")))
    year = str(date["date-parts"][0][0])
    title = " ".join(metadata.get("title", []))
    if metadata.get("subtitle"):
        title += ": " + " ".join(metadata["subtitle"])
    return {"key": entry.key, "doi": doi, "title": field(entry, "title"), "source_title": title,
            "title_match": normalize(field(entry, "title")) == normalize(title),
            "year": field(entry, "year"), "source_year": year,
            "volume": field(entry, "volume"), "source_volume": metadata.get("volume", ""),
            "number": field(entry, "number"), "source_number": metadata.get("issue", ""),
            "pages": field(entry, "pages"), "source_pages": metadata.get("page", metadata.get("article-number", "")),
            "source_authors": "; ".join(" ".join([author.get("given", ""), author.get("family", "")]).strip()
                                        for author in metadata.get("author", [])),
            "source_url": metadata.get("URL", ""), "metadata_file": str(path.relative_to(root))}


entries = [entry for entry in bibliography.entries if field(entry, "doi")]
records = []
with ThreadPoolExecutor(max_workers=4) as executor:
    futures = {executor.submit(fetch, entry): entry.key for entry in entries}
    for future in as_completed(futures):
        records.append(future.result())
        print("Verified DOI", futures[future], flush=True)
records.sort(key=lambda record: record["key"])
with (destination / "comparison.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(records[0]))
    writer.writeheader()
    writer.writerows(records)
for record in records:
    differences = [name for name in ("year", "volume", "number", "pages")
                   if record[f"source_{name}"] and normalize(record[name]) != normalize(record[f"source_{name}"])]
    if differences or not record["title_match"]:
        print(json.dumps({"key": record["key"], "title_match": record["title_match"],
                          "differences": {name: [record[name], record[f"source_{name}"]] for name in differences}},
                         ensure_ascii=False))
print("Verified", len(records), "DOI records")
