from pathlib import Path
from pypdf import PdfReader

paper = Path("../research/jmlr221217/22-1217.pdf")
reader = PdfReader(paper)
output = Path("results/scratch/submission_sqr_source.txt")
output.write_text("\n\n".join(f"PAGE {index + 1}\n{page.extract_text()}" for index, page in enumerate(reader.pages)))
print("Read", len(reader.pages), "pages")
