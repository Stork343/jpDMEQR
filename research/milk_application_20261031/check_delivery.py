import csv
import json
from datetime import datetime, timezone
from pathlib import Path

import pymupdf


root = Path.cwd()
output = root / "results/preflight/milk_scs_20261031"
completion = json.loads((output / "completion_manifest.json").read_text())
with (output / "inference.csv").open(newline="") as stream:
    inference = list(csv.DictReader(stream))
interpretation = output / "interpretation"
with (interpretation / "quantile_contrasts.csv").open(newline="") as stream:
    contrasts = list(csv.DictReader(stream))
with (interpretation / "quantile_equality_test.csv").open(newline="") as stream:
    global_test = list(csv.DictReader(stream))[0]
with (interpretation / "observed_distribution.csv").open(newline="") as stream:
    descriptive = list(csv.DictReader(stream))
if len(inference) != completion["cases_registered"] or len(inference) != 18:
    raise ValueError("The empirical result matrix is incomplete.")
run = completion["execution_version"]
backed_up = 0
for row in inference:
    source = root / "results/raw/application/milk_scs_20261031" / run / row["case_id"] / "fit.rds"
    backup = Path("/Users/houjian/jpDMEQR-backups/milk_scs_20261031") / run / row["case_id"] / "fit.rds"
    if source.stat().st_size == 0 or source.stat().st_size != backup.stat().st_size:
        raise ValueError("An empirical fit is absent from its independent backup.")
    backed_up += 1
documents = []
for name, heading, rows in (
    ("sjs_profile_quantile_main", "Milk spectra and somatic cell score", [row for row in inference if row["setting"] == "reference"]),
    ("sjs_profile_quantile_supplement", "Milk spectra construction and analysis", inference),
):
    path = root / "manuscript" / f"{name}.pdf"
    document = pymupdf.open(path)
    text = " ".join(" ".join(page.get_text() for page in document).split())
    text = text.translate(str.maketrans({"−": "-", "–": "-"}))
    normalised = text.casefold()
    if heading.casefold() not in normalised or any(term in normalised for term in ("gse65391", "plasma-cell", "??")):
        raise ValueError("The compiled application text does not match the current study.")
    application = normalised[normalised.index(heading.casefold()):]
    for row in rows:
        for field in ("beta_tilde", "simultaneous_lower", "simultaneous_upper"):
            value = f"{float(row[field]):.3f}"
            if value not in application:
                raise ValueError(f"A saved empirical value is absent from {name}: {value}")
    if f"{float(global_test['p_value']):.3f}" not in application:
        raise ValueError("The quantile comparison p value is absent from the compiled text.")
    extra_checks = 0
    if name == "sjs_profile_quantile_supplement":
        for row in contrasts:
            for field in ("estimate", "se", "simultaneous_lower", "simultaneous_upper", "bonferroni_p"):
                if f"{float(row[field]):.3f}" not in application:
                    raise ValueError("A contrast value is absent from the compiled supplement.")
                extra_checks += 1
        caption = "observed somatic cell score distributions"
        if caption not in normalised:
            raise ValueError("The observed-distribution table caption is absent.")
        descriptive_text = normalised[normalised.index(caption):]
        group_labels = {"all": "all eligible cows", "A1": "a1 milk type", "A2": "a2 milk type",
                        "single_control": "one milk control", "two_controls": "two milk controls"}
        for row in descriptive:
            values = [group_labels[row["group"]], row["cows"], row["samples"]]
            values.extend(f"{float(row[field]):.3f}" for field in
                          ("cow_weighted_mean", "cow_weighted_q25", "cow_weighted_q50", "cow_weighted_q75"))
            if " ".join(values) not in descriptive_text:
                raise ValueError("A complete observed-distribution row is absent from the compiled supplement.")
            extra_checks += 4
    overflow = []
    for page_number, page in enumerate(document, start=1):
        for word in page.get_text("words"):
            if word[0] < -0.5 or word[1] < -0.5 or word[2] > page.rect.width + 0.5 or word[3] > page.rect.height + 0.5:
                overflow.append(page_number)
    if overflow:
        raise ValueError("Extracted PDF text extends beyond its page boundary.")
    documents.append({"file": str(path.relative_to(root)), "pages": len(document),
                      "verified_result_rows": len(rows), "additional_values_checked": extra_checks,
                      "out_of_page_text": 0})
report = {"verified_utc": datetime.now(timezone.utc).isoformat(), "completed_cases": len(inference),
          "backed_up_fits": backed_up, "documents": documents, "review_method": "PDF text extraction and page coordinates"}
(output / "delivery_checks.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
print(json.dumps(report, ensure_ascii=False, indent=2))
