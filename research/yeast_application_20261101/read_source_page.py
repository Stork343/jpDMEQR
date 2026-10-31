from pathlib import Path

from bs4 import BeautifulSoup

root = Path("results/scratch/empirical_selection_20261101")
path = root / "yeast_2019_article.html"
soup = BeautifulSoup(path.read_text(encoding="utf-8"), "html.parser")
for node in soup(["script", "style", "nav"]):
    node.decompose()
text = soup.get_text("\n", strip=True)
(root / "yeast_2019_article.txt").write_text(text, encoding="utf-8")
for paragraph in soup.find_all(["p", "h2", "h3"]):
    value = paragraph.get_text(" ", strip=True)
    if any(term in value.lower() for term in ("tetrad", "4390", "4,390", "by-rm", "flow cyt", "random spor", "dissect", "data availability", "duplicate")):
        print(value)
