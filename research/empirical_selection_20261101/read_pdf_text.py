import argparse
from pathlib import Path

from pypdf import PdfReader

parser = argparse.ArgumentParser()
parser.add_argument("path", type=Path)
args = parser.parse_args()
reader = PdfReader(args.path)
destination = args.path.with_suffix(".txt")
destination.write_text("\n".join(page.extract_text() for page in reader.pages), encoding="utf-8")
print(len(reader.pages), destination)
