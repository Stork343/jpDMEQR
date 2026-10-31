import csv
import json
from pathlib import Path

import pymupdf


root = Path(__file__).resolve().parents[2]
tables = root / "results/tables/soynam_20261031"
manuscript = root / "manuscript"


def normalise(text):
    return " ".join(text.replace("−", "-").replace("–", "-").split())


def read_rows(name):
    with (tables / name).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


records = []
texts = {}
for name in ("sjs_profile_quantile_main", "sjs_profile_quantile_supplement"):
    path = manuscript / f"{name}.pdf"
    document = pymupdf.open(path)
    pages = []
    outside = []
    for page_index, page in enumerate(document):
        pages.append(normalise(page.get_text()))
        for word in page.get_text("words"):
            if word[0] < -1 or word[1] < -1 or word[2] > page.rect.width + 1 or word[3] > page.rect.height + 1:
                outside.append({"page": page_index + 1, "word": word[4]})
    text = " ".join(pages)
    if outside:
        raise ValueError(f"Text extends beyond the PDF page: {name}: {outside[:5]}")
    if "??" in text:
        raise ValueError(f"An unresolved reference remains in {name}")
    if "soybean" not in text.lower() or "5533" not in text.replace(",", ""):
        raise ValueError(f"The updated application is absent from {name}")
    if "somatic cell score" in text.lower():
        raise ValueError(f"The application text is inconsistent in {name}")
    texts[name] = text
    records.append({"file": str(path), "pages": len(document), "characters": len(text), "out_of_page_words": len(outside)})

supplement = texts["sjs_profile_quantile_supplement"]
for row in read_rows("inference.csv"):
    for column in ("estimate_bu_ac", "simultaneous_lower", "simultaneous_upper"):
        value = f"{float(row[column]):.3f}"
        if value not in supplement:
            raise ValueError(f"A sensitivity value is absent from the supplement: {row['case_id']} {column} {value}")

for row in read_rows("reference_contrasts.csv"):
    for column in ("estimate", "se", "lower", "upper"):
        value = f"{float(row[column]):.3f}"
        if value not in supplement:
            raise ValueError(f"A coordinate or contrast value is absent: {row['parameter']} {column} {value}")

for row in read_rows("prediction_metrics.csv"):
    value = f"{float(row['loss']):.3f}"
    if value not in texts["sjs_profile_quantile_main"]:
        raise ValueError(f"A prediction loss is absent from the article: {row['method']} {row['tau']} {value}")

figure = pymupdf.open(manuscript / "generated/soynam_profile/results.pdf")
figure_text = normalise(figure[0].get_text())
for phrase in ("Canopy association", "Held-out prediction", "Quantile", "2013", "2014"):
    if phrase not in figure_text:
        raise ValueError(f"A figure label is absent: {phrase}")

result = {"status": "passed", "documents": records, "inference_rows": 42,
          "joint_parameters": 11, "prediction_rows": 9,
          "inspection": "PDF text and page coordinates"}
destination = root / "results/preflight/soynam_20261031/pdf_checks.json"
destination.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(json.dumps(result, ensure_ascii=False, indent=2))
