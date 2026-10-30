import csv
import gzip
from pathlib import Path


for release in ("R1", "R2"):
    path = Path("data/raw/GSE65391") / f"GSE65391_non-normalized_data_Illumina_HT12_V4_{release}.txt.gz.partial"
    with gzip.open(path, "rt") as stream:
        reader = csv.reader(stream, delimiter="\t")
        header = next(reader)
        first = next(reader)
    print(release, "columns", len(header), "row_columns", len(first))
    print("header_begin", header[:12])
    print("header_end", header[-6:])
    print("first_identifier", first[0])
