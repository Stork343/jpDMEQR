import json
from datetime import datetime, timezone
from pathlib import Path
from zipfile import ZipFile


metadata_path = Path("results/scratch/dataset_search_20261030/milk_metadata.json")
metadata = json.loads(metadata_path.read_text())
version = metadata["data"]["latestVersion"]
local_files = {
    296867: Path("data/raw/candidate_datasets_20261030/Database_MIR_A2_milk.xlsx"),
    320914: Path("data/raw/candidate_datasets_20261030/milk_README.txt"),
}
files = []
for file_id, local_path in local_files.items():
    matches = [entry for entry in version["files"] if entry["dataFile"]["id"] == file_id]
    if len(matches) != 1:
        raise ValueError("The deposited file identifier is not unique.")
    entry = matches[0]
    if entry["restricted"] or local_path.stat().st_size != entry["dataFile"]["filesize"]:
        raise ValueError("The source file does not match the public repository record.")
    if local_path.suffix == ".xlsx":
        with ZipFile(local_path) as archive:
            if archive.testzip() is not None:
                raise ValueError("The Excel archive failed its integrity check.")
    files.append({"file_id": file_id, "filename": entry["dataFile"]["filename"],
                  "local_path": str(local_path), "bytes": local_path.stat().st_size,
                  "url": f"https://dataverse.csuc.cat/api/access/datafile/{file_id}"})
report = {"dataset_doi": "10.34810/data2338",
          "metadata_url": "https://dataverse.csuc.cat/api/datasets/:persistentId?persistentId=doi:10.34810/data2338",
          "version": f"{version['versionNumber']}.{version['versionMinorNumber']}",
          "release_time": version["releaseTime"], "license": version["license"]["name"],
          "files": files, "verified_utc": datetime.now(timezone.utc).isoformat()}
output = Path("results/preflight/dataset_search_20261030/milk_source_manifest.json")
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
print(json.dumps(report, ensure_ascii=False, indent=2))
