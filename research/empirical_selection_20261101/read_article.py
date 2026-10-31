import argparse
from pathlib import Path

from bs4 import BeautifulSoup

parser = argparse.ArgumentParser()
parser.add_argument("path", type=Path)
args = parser.parse_args()
soup = BeautifulSoup(args.path.read_text(encoding="utf-8"), "html.parser")
article = soup.find("article") or soup
for node in article(["script", "style", "nav"]):
    node.decompose()
text = article.get_text("\n", strip=True)
destination = args.path.with_suffix(".txt")
destination.write_text(text, encoding="utf-8")
print(destination)
