import argparse
from pathlib import Path
from urllib.parse import urljoin

from bs4 import BeautifulSoup

parser = argparse.ArgumentParser()
parser.add_argument("path", type=Path)
parser.add_argument("base_url")
args = parser.parse_args()
soup = BeautifulSoup(args.path.read_text(encoding="utf-8"), "html.parser")
for link in soup.find_all("a", href=True):
    if "MediaObjects" in link["href"] or "supplement" in link["href"].lower():
        print(link.get_text(" ", strip=True), urljoin(args.base_url, link["href"]))
