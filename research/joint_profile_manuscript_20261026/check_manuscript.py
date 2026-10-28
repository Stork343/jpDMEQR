import json
from pathlib import Path
import re
import subprocess

import bibtexparser
import pdfplumber
from pypdf import PdfReader
from pylatexenc.latexwalker import LatexWalker, LatexMacroNode, LatexGroupNode, LatexEnvironmentNode, LatexMathNode, get_default_latex_context_db
from pylatexenc.macrospec import MacroSpec

root = Path.cwd()
manuscript = root / "manuscript"
output = root / "results/scratch/joint_profile_manuscript_20261026/pdf_checks"
output.mkdir(parents=True, exist_ok=True)
context = get_default_latex_context_db()
context.add_context_category("manuscript", macros=[MacroSpec(name, "{") for name in
    ["input", "include", "label", "ref", "eqref", "citeA"]] + [MacroSpec("cite", "[{")], prepend=True)
bibliography = bibtexparser.parse_file(manuscript / "references.bib")
if bibliography.failed_blocks:
    raise ValueError("Unparsed bibliography entries")
bibkeys = {entry.key for entry in bibliography.entries}


def nodes_walk(nodes):
    for node in nodes:
        yield node
        if isinstance(node, (LatexGroupNode, LatexEnvironmentNode, LatexMathNode)):
            yield from nodes_walk(node.nodelist)
        if isinstance(node, LatexMacroNode) and node.nodeargd:
            for argument in node.nodeargd.argnlist:
                if argument is not None:
                    yield from nodes_walk([argument])


def source_check(entry):
    citations = set()
    labels = set()
    references = set()
    expanded = subprocess.run(["/Library/TeX/texbin/latexpand", "--fatal", entry.name],
                              cwd=manuscript, check=True, capture_output=True, text=True).stdout
    nodes, _, _ = LatexWalker(expanded, latex_context=context, tolerant_parsing=False).get_latex_nodes()
    for node in nodes_walk(nodes):
        if not isinstance(node, LatexMacroNode) or not node.nodeargd:
            continue
        arguments = [arg for arg in node.nodeargd.argnlist if isinstance(arg, LatexGroupNode)]
        if not arguments:
            continue
        value = arguments[-1].latex_verbatim()[1:-1]
        if node.macroname in ("cite", "citeA"):
            citations.update(key.strip() for key in value.split(","))
        elif node.macroname == "label":
            if value in labels:
                raise ValueError(f"Repeated label {value}")
            labels.add(value)
        elif node.macroname in ("ref", "eqref"):
            references.add(value)
    missing_citations = sorted(citations - bibkeys)
    missing_references = sorted(references - labels)
    if missing_citations or missing_references:
        raise ValueError(f"Unresolved source links {missing_citations} {missing_references}")
    return {"expanded_characters": len(expanded), "citations": len(citations), "labels": len(labels)}


def font_check(font):
    font = font.get_object()
    if "/DescendantFonts" in font:
        return all(font_check(item) for item in font["/DescendantFonts"])
    if font.get("/Subtype") == "/Type3":
        return True
    descriptor = font.get("/FontDescriptor")
    if descriptor is None:
        return False
    return any(key in descriptor.get_object() for key in ("/FontFile", "/FontFile2", "/FontFile3"))


report = {"inspection": "PDF object and text inspection without visual inspection", "documents": []}
for name in ("sjs_profile_quantile_main", "sjs_profile_quantile_supplement"):
    source = source_check(manuscript / f"{name}.tex")
    reader = PdfReader(manuscript / f"{name}.pdf")
    fonts = {}
    for page in reader.pages:
        resources = page["/Resources"].get_object()
        for font in resources["/Font"].get_object().values():
            value = font.get_object()
            fonts[str(value.get("/BaseFont"))] = font_check(font)
    if not all(fonts.values()):
        raise ValueError(f"Unembedded fonts in {name}: {fonts}")
    pages = []
    page_text = []
    with pdfplumber.open(manuscript / f"{name}.pdf") as document:
        for number, page in enumerate(document.pages, 1):
            chars = [char for char in page.chars if char["text"].strip()]
            outside = [char for char in chars if char["x0"] < -1 or char["x1"] > page.width + 1
                       or char["top"] < -1 or char["bottom"] > page.height + 1]
            text = page.extract_text() or ""
            if not chars or "??" in text:
                raise ValueError(f"Empty or unresolved PDF page {name} {number}")
            pages.append({"page": number, "characters": len(chars), "outside": outside,
                          "x_min": min(char["x0"] for char in chars), "x_max": max(char["x1"] for char in chars),
                          "y_min": min(char["top"] for char in chars), "y_max": max(char["bottom"] for char in chars)})
            page_text.append(text)
    (output / f"{name}.txt").write_text("\n\n".join(page_text))
    log = (manuscript / f"{name}.log").read_text()
    bad_log = [line for line in log.splitlines() if re.search(r"undefined|multiply defined|Missing character|^!", line)]
    if bad_log:
        raise ValueError(f"LaTeX diagnostics in {name}: {bad_log}")
    report["documents"].append({"name": name, "pages": len(pages), "source": source,
                                "fonts": fonts, "page_bounds": pages,
                                "title": str(reader.metadata.title)})
(output / "report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False))
for item in report["documents"]:
    outside_count = sum(len(page["outside"]) for page in item["page_bounds"])
    print(item["name"], "pages", item["pages"], "outside characters", outside_count,
          "embedded fonts", len(item["fonts"]), "source", item["source"])
    for page in item["page_bounds"]:
        if page["outside"]:
            print("Page", page["page"], "outside text", "".join(char["text"] for char in page["outside"]))
if any(page["outside"] for item in report["documents"] for page in item["page_bounds"]):
    raise ValueError("PDF text extends outside its page")
