import json
import zipfile
from pathlib import Path

import pandas as pd

root = Path("data/raw/covid_charite_20261105")
metadata = pd.read_excel(root / "Table_S9.xlsx")
measurements = pd.read_excel(root / "Table_S3.xlsx")
assert not metadata.File.duplicated().any()
assert not metadata[["Patient Id", "Time point - days since admission"]].duplicated().any()
assert metadata["Patient Id"].notna().all()
available = metadata.dropna(subset=["Horowitz index"])
assert available["Horowitz index"].gt(0).all()
matrix_path = root / "Table_S8_complete.xlsx"
matrix_available = matrix_path.exists() and matrix_path.stat().st_size == 2735423 and zipfile.is_zipfile(matrix_path)
if matrix_available:
    with zipfile.ZipFile(matrix_path) as archive:
        assert archive.testzip() is None
output = Path("results/preflight/empirical_selection_20261101/covid_charite_metadata")
output.mkdir(parents=True, exist_ok=True)
summary = {
    "phase": "metadata_audit",
    "source_doi": "10.1016/j.cels.2021.05.005",
    "participants": int(metadata["Patient Id"].nunique()),
    "sampling_points": len(metadata),
    "maximum_sampling_points_per_participant": int(metadata.groupby("Patient Id").size().max()),
    "measurement_definitions": len(measurements),
    "horowitz_records": len(available),
    "horowitz_participants": int(available["Patient Id"].nunique()),
    "horowitz_minimum": float(available["Horowitz index"].min()),
    "horowitz_maximum": float(available["Horowitz index"].max()),
    "horowitz_above_1000": int(available["Horowitz index"].gt(1000).sum()),
    "bmi_missing_among_horowitz_records": int(available.BMI.isna().sum()),
    "source_clinical_timing": "latest_value_from_up_to_two_preceding_days_when_current_day_missing",
    "individual_clinical_measurement_timestamps_provided": False,
    "protein_matrix_complete": bool(matrix_available),
    "roles_assigned": False,
    "outcome_associations_fitted": False,
    "formal_validation_used": False,
    "manuscript_adoption": False,
}
(output / "status.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
metadata.groupby("Patient Id").size().value_counts().sort_index().rename_axis("sampling_points").rename("participants").to_csv(output / "repeat_counts.csv")
print(json.dumps(summary, indent=2))
