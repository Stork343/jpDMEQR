from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
from dulwich.lfs import LFSPointer

root = Path("data/raw/mpc_20261106")
root.mkdir(parents=True, exist_ok=True)
commit = "a8ea891da32add81fabaa8e744cf8cb72a8a0e66"
clinical = {"mpcproject_broad_2021": ["data_clinical_patient.txt", "data_clinical_sample.txt", "data_timeline_psa.txt",
    "data_timeline_treatment.txt", "data_timeline_sample.txt", "data_gene_panel_matrix.txt"],
    "msk_chord_2024": ["data_clinical_patient.txt", "data_clinical_sample.txt", "data_timeline_psa_labs.txt",
    "data_timeline_treatment.txt", "data_timeline_specimen.txt", "data_gene_panel_matrix.txt", "data_timeline_gleason.txt",
    "data_timeline_cea_labs.txt","data_timeline_ca_19-9_labs.txt","data_timeline_ca_15-3_labs.txt"]}


def fetch(url, dest):
    partial = dest.with_name(dest.name + ".incoming")
    subprocess.run(["curl", "--http1.1", "--silent", "--show-error", "--fail", "--location", "--retry", "5",
                    "--retry-all-errors", "--connect-timeout", "12", "--max-time", "75", "--proxy", "http://127.0.0.1:7897",
                    "--output", str(partial), url], check=True)
    partial.replace(dest)


def obtain(job):
    study, filename = job
    dest = root / study
    dest.mkdir(exist_ok=True)
    pointers = dest / "published_pointers"
    pointers.mkdir(exist_ok=True)
    pointer = pointers / filename
    suffix = f"cBioPortal/datahub/{commit}/public/{study}/{filename}"
    if not pointer.exists():
        fetch("https://raw.githubusercontent.com/" + suffix, pointer)
    fields = LFSPointer.from_bytes(pointer.read_bytes())
    if fields is None:
        raise ValueError(f"无效的 Git LFS 来源 {pointer}")
    size, digest = fields.size, fields.oid
    record = {"study": study, "file": filename, "bytes": size, "sha256": digest}
    if mode == "index":
        return record
    path = dest / filename
    if not path.exists():
        with tarfile.open(root / f"{study}.tar.gz","r:gz") as archive:
            member, = [m for m in archive.getmembers() if m.isfile() and Path(m.name).name == filename]
            assert member.size == size
            partial = path.with_name(path.name + ".archive.incoming")
            with archive.extractfile(member) as stream, partial.open("wb") as output:
                shutil.copyfileobj(stream, output)
            with partial.open("rb") as stream:
                assert hashlib.file_digest(stream,"sha256").hexdigest() == digest
            partial.replace(path)
    else:
        with path.open("rb") as stream:
            assert path.stat().st_size == size and hashlib.file_digest(stream,"sha256").hexdigest() == digest
    print(study, filename, size, flush=True)
    return record


mode = sys.argv[1]
assert mode in {"index","clinical","genomics"}
archives = {"msk_chord_2024": 46876803, "mpcproject_broad_2021": 22961713}
if mode != "index":
    for study, size in archives.items():
        path = root / f"{study}.tar.gz"
        if not path.with_name(path.name + ".verified.json").exists():
            subprocess.run([sys.executable,"research/loop_application_20261106/download.py",
                f"https://datahub.assets.cbioportal.org/{study}.tar.gz",str(path),str(size),"gzip"],
                env=dict(os.environ, CLINICAL_DOWNLOAD_WORKERS="8"),check=True)
jobs = [(s,n) for s, files in clinical.items() for n in (files if mode == "clinical" else
        [*files,"data_cna.txt","data_mutations.txt"] if mode == "index" else ["data_cna.txt","data_mutations.txt"])]
with ThreadPoolExecutor(max_workers=3) as executor:
    results = list(executor.map(obtain, jobs))
(root / f"{mode}_manifest.json").write_text(json.dumps(results,indent=2)+"\n")
