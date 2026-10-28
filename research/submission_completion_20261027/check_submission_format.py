import json
from pathlib import Path
import subprocess

from pylatexenc.latex2text import LatexNodes2Text
from pylatexenc.latexwalker import LatexWalker, LatexMacroNode, LatexGroupNode, LatexMathNode, LatexEnvironmentNode, get_default_latex_context_db
from pylatexenc.macrospec import MacroSpec

root = Path.cwd()
source = subprocess.run(["/Library/TeX/texbin/latexpand", "--fatal", "sjs_profile_quantile_main.tex"],
                        cwd=root / "manuscript", check=True, capture_output=True, text=True).stdout
context = get_default_latex_context_db()
context.add_context_category("submission", macros=[MacroSpec("abstract", "[{"), MacroSpec("keywords", "{"),
                            MacroSpec("titlemark", "{")], prepend=True)
nodes, _, _ = LatexWalker(source, latex_context=context, tolerant_parsing=False).get_latex_nodes()
values = {}


def walk(items):
    for item in items:
        yield item
        if isinstance(item, (LatexEnvironmentNode, LatexGroupNode, LatexMathNode)):
            yield from walk(item.nodelist)


for node in walk(nodes):
    if isinstance(node, LatexMacroNode) and node.macroname in {"abstract", "keywords", "titlemark"}:
        arguments = [item for item in node.nodeargd.argnlist if isinstance(item, LatexGroupNode)]
        values[node.macroname] = arguments[-1]
if set(values) != {"abstract", "keywords", "titlemark"}:
    raise ValueError("The article submission fields are incomplete")
abstract = values["abstract"]
if any(isinstance(node, LatexMathNode) for node in abstract.nodelist):
    raise ValueError("The abstract contains a mathematical formula")
words = len(LatexNodes2Text().nodelist_to_text(abstract.nodelist).split())
numeric_macros = {"JointIntervals", "JointDatasets", "JointRuntimeTenThousand", "RegisteredJointDatasets", "RegisteredJointSettings"}
words += sum(isinstance(node, LatexMacroNode) and node.macroname in numeric_macros for node in abstract.nodelist)
keywords = [text.strip() for text in LatexNodes2Text().nodelist_to_text(values["keywords"].nodelist).split(";")]
running_title = LatexNodes2Text().nodelist_to_text(values["titlemark"].nodelist)
if words > 150 or len(keywords) > 7 or len(running_title) > 40:
    raise ValueError(f"A journal field exceeds its limit: {words}, {len(keywords)}, {len(running_title)}")
if keywords != sorted(keywords, key=str.casefold):
    raise ValueError("The keywords are not alphabetically ordered")
record = {"abstract_words": words, "keywords": keywords, "running_title": running_title,
          "running_title_characters": len(running_title),
          "source": "https://onlinelibrary.wiley.com/page/journal/14679469/homepage/forauthors.html",
          "limits": {"abstract_words": 150, "keywords": 7, "running_title_characters": 40}}
destination = root / "results/scratch/submission_completion_20261027/submission_format.json"
destination.write_text(json.dumps(record, ensure_ascii=False, indent=2))
print(json.dumps(record, ensure_ascii=False, indent=2))
