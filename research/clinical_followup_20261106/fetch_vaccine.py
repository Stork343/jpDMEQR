from pathlib import Path
from time import time
import shutil

from remotezip import RemoteZip

folder = Path("data/raw/vaccine_20261106")
folder.mkdir(exist_ok=True, parents=True)
url = "https://ndownloader.figshare.com/files/31615370?retrieval="+str(int(time()))
with RemoteZip(url, timeout=180) as archive:
    for filename in ["dataset_metadata.csv", "immdata_all.rds", "all_noNorm_withResponse_eset.rds"]:
        member = archive.getinfo("data/"+filename)
        destination = folder/filename
        if not destination.exists():
            partial = destination.with_suffix(destination.suffix+".incoming")
            with archive.open(member) as origin, partial.open("wb") as stream:
                shutil.copyfileobj(origin, stream, 1024*1024)
            assert partial.stat().st_size == member.file_size
            partial.replace(destination)
        assert destination.stat().st_size == member.file_size
        print(filename, member.file_size, flush=True)
