import json
from pathlib import Path

import pandas as pd

root = Path("data/raw/drugmatrix_20261101")
output = Path("results/preflight/empirical_selection_20261101/drugmatrix_metadata")
output.mkdir(parents=True, exist_ok=True)
clinical = pd.read_excel(root / "DrugMatrix_ClinicalChemistry_Hematology.xlsx", sheet_name="DM_ClinicalChemistry_Hematology")
clinical.columns = clinical.columns.str.strip()
labels = pd.read_csv(root / "Codelink_Microarray_Data_Labels_spreadsheet_02_23_2012.txt", sep="\t", low_memory=False)
labels.columns = labels.columns.str.strip()
assert not clinical["ANIMAL_ID"].duplicated().any()
assert not labels["ARRAY_ID"].duplicated().any()
liver = labels.loc[labels["ORGAN_OR_CELL_TYPE"].str.upper().eq("LIVER") & labels["STUDY_TYPE"].eq("In Vivo")].copy()
matched = liver.merge(clinical, on="ANIMAL_ID", how="left", validate="many_to_one", suffixes=("_array", "_clinical"), indicator=True)
treated = matched.loc[matched["TEST_ARTICLE"].notna()].copy()
alt = pd.to_numeric(treated["ALT"], errors="coerce")
invalid = treated.loc[treated["ALT"].notna() & alt.isna(), "ALT"].value_counts()
invalid.to_csv(output / "nonnumeric_alt_codes.csv")
eligible = treated.loc[alt.notna() & alt.gt(0)].copy()
eligible["ALT_numeric"] = alt.loc[eligible.index]
eligible = eligible.sort_values("ARRAY_ID").drop_duplicates("ANIMAL_ID")
compounds = eligible.groupby("TEST_ARTICLE").agg(animals=("ANIMAL_ID", "nunique"), studies=("STUDY_ID", "nunique"))
studies = eligible.groupby("STUDY_ID").agg(animals=("ANIMAL_ID", "nunique"), compounds=("TEST_ARTICLE", "nunique"))
compounds.to_csv(output / "compound_structure.csv")
studies.to_csv(output / "study_structure.csv")
summary = {
    "clinical_animals": len(clinical),
    "codelink_arrays": len(labels),
    "liver_in_vivo_arrays": len(liver),
    "arrays_matched_to_clinical_animals": int(matched["_merge"].eq("both").sum()),
    "treated_arrays_with_positive_alt": int((alt.notna() & alt.gt(0)).sum()),
    "eligible_unique_animals": len(eligible),
    "compounds": len(compounds),
    "studies": len(studies),
    "maximum_animals_per_compound": int(compounds["animals"].max()),
    "maximum_animals_per_study": int(studies["animals"].max()),
    "compounds_observed_in_multiple_studies": int(compounds["studies"].gt(1).sum()),
    "studies_with_multiple_compounds": int(studies["compounds"].gt(1).sum()),
    "distinct_alt_values": int(eligible["ALT_numeric"].nunique()),
    "outcome_associations_fitted": False,
    "manuscript_adoption": False,
}
(output / "audit.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
eligible.to_csv(root / "matched_liver_animals.csv", index=False)
print(json.dumps(summary, indent=2))
print(clinical["TREATMENT"].value_counts().to_string())
print(eligible["ALT_numeric"].describe().to_string())
