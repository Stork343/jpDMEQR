import json
from pathlib import Path

import biom
import numpy as np
import pandas as pd

raw = Path("data/raw/ibdmdb_20261101")
output = Path("data/derived/ibdmdb_20261105")
report = Path("results/preflight/empirical_selection_20261101/ibdmdb_cohort")
output.mkdir(parents=True, exist_ok=True)
report.mkdir(parents=True, exist_ok=True)
metadata = pd.read_csv(raw / "hmp2_metadata_2018-08-20.csv", low_memory=False)
measurements = pd.concat([
    metadata.loc[metadata[column].notna(), ["site_sub_coll", "Participant ID", column]]
    .rename(columns={column: "value"}).assign(source_field=column)
    for column in ["fecalcal", "fecalcal_ng_ml"]
], ignore_index=True)
assert measurements.value.gt(0).all()
identity = measurements.groupby("site_sub_coll")["Participant ID"].nunique()
assert identity.eq(1).all()
agreement = measurements.groupby("site_sub_coll").agg(
    minimum=("value", "min"), maximum=("value", "max"), records=("value", "size"),
    fields=("source_field", "nunique"),
)
agreement["range"] = agreement.maximum - agreement.minimum
assert agreement["range"].max() < 0.001
measurements["priority"] = measurements.source_field.map({"fecalcal_ng_ml": 0, "fecalcal": 1})
canonical = measurements.sort_values(["site_sub_coll", "priority"], kind="stable").drop_duplicates("site_sub_coll")
canonical = canonical.rename(columns={"value": "reported_calprotectin"})
metabolomics = metadata.loc[metadata.data_type.eq("metabolomics")].copy()
assert not metabolomics["External ID"].duplicated().any()
assert not metabolomics.site_sub_coll.duplicated().any()
cohort = metabolomics.merge(canonical[["site_sub_coll", "reported_calprotectin", "source_field"]],
                            on="site_sub_coll", how="inner", validate="one_to_one")
required = ["Participant ID", "External ID", "site_sub_coll", "consent_age", "sex", "diagnosis", "site_name", "week_num",
            "Antibiotics", "Immunosuppressants (e.g. oral corticosteroids)"]
assert not cohort[required].isna().any().any()
assert cohort.consent_age.between(0, 100).all()
assert set(cohort.Antibiotics).issubset({"Yes", "No"})
assert set(cohort["Immunosuppressants (e.g. oral corticosteroids)"]).issubset({"Yes", "No"})
assert set(cohort.sex).issubset({"Female", "Male"})
assert set(cohort.diagnosis).issubset({"CD", "UC", "nonIBD"})
participants = cohort[["Participant ID", "diagnosis", "site_name", "sex", "consent_age"]].drop_duplicates()
assert not participants["Participant ID"].duplicated().any()
participants = participants.sort_values("Participant ID").reset_index(drop=True)
rng = np.random.default_rng(20261105)
participants["role"] = "discovery"
participants["outer_fold"] = 0
for site, indices in participants.groupby("site_name", sort=True).groups.items():
    order = rng.permutation(np.asarray(indices))
    reserved = max(1, round(0.2 * len(order)))
    participants.loc[order[:reserved], "role"] = "validation_reserved"
    discovery = order[reserved:]
    participants.loc[discovery, "outer_fold"] = np.arange(len(discovery)) % 5 + 1
if (output / "roles.csv").exists():
    assert pd.read_csv(output / "roles.csv").equals(participants)
else:
    participants.to_csv(output / "roles.csv", index=False)
cohort = cohort.merge(participants[["Participant ID", "role", "outer_fold"]], on="Participant ID", validate="many_to_one")
rows = cohort[required + ["reported_calprotectin", "source_field", "role", "outer_fold"]].rename(columns={
    "Participant ID": "strain_id", "External ID": "sample_id", "consent_age": "age", "site_name": "site",
    "Antibiotics": "antibiotics", "Immunosuppressants (e.g. oral corticosteroids)": "immunosuppressants",
})
rows["y"] = np.log(rows.reported_calprotectin)
rows["condition"] = "fecalcal"
rows = rows.sort_values(["strain_id", "week_num", "sample_id"])
for role in ["discovery", "validation_reserved"]:
    rows.loc[rows.role.eq(role)].to_csv(output / f"{role}.csv", index=False)

table = biom.load_table(str(raw / "HMP2_metabolomics_w_metadata.biom"))
annotation = pd.read_csv(raw / "feature_annotation.csv")
named = annotation.loc[annotation.Metabolite.notna()].copy()
assert named["Pooled QC sample CV"].notna().all()
assert named.Method.notna().all()
named = named.sort_values(["Metabolite", "Pooled QC sample CV", "feature_id"], kind="stable").drop_duplicates("Metabolite")
named = named.loc[named["Pooled QC sample CV"].le(0.3)].sort_values("feature_id")
named.to_csv(output / "feature_annotation.csv", index=False)
sample_ids = list(table.ids())
assert set(rows.sample_id).issubset(sample_ids)
feature_lookup = pd.Series(np.arange(table.shape[0]), index=table.ids(axis="observation"))
values = table.matrix_data[feature_lookup.loc[named.feature_id].to_numpy(), :].toarray().T
assert np.isfinite(values).all() and (values >= 0).all()
pd.DataFrame(values, index=sample_ids, columns=named.feature_id).to_csv(
    output / "metabolite_abundance.csv.gz", index_label="sample_id")
signal = pd.DataFrame(index=sample_ids)
for method in sorted(annotation.Method.dropna().unique()):
    current = annotation.loc[annotation.Method.eq(method), "feature_id"]
    matrix = table.matrix_data[feature_lookup.loc[current].to_numpy(), :].toarray().T
    medians = np.array([np.median(row[row > 0]) for row in matrix])
    assert np.isfinite(medians).all() and (medians > 0).all()
    signal[method] = np.log2(medians)
signal.to_csv(output / "method_signal.csv", index_label="sample_id")
rows.groupby(["role", "diagnosis"]).agg(participants=("strain_id", "nunique"), records=("sample_id", "size")).to_csv(report / "cohort_counts.csv")
participants.groupby(["role", "site_name"]).size().rename("participants").to_csv(report / "site_counts.csv")
summary = {
    "status": "prepared", "participants": len(participants), "records": len(rows),
    "named_features": len(named), "maximum_measurements": int(rows.groupby("strain_id").size().max()),
    "features_without_method_annotation": int(annotation.Method.isna().sum()),
    "exact_specimen_matching": True, "matching_key": "site_sub_coll",
    "specimens_with_both_source_fields": int(agreement.fields.eq(2).sum()),
    "maximum_source_field_difference": float(agreement["range"].max()),
    "source_units": "ng_per_ml_in_fecalcal_ng_ml_export",
    "outcome_associations_fitted": False, "formal_validation_used": False, "manuscript_adoption": False,
}
(report / "status.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
print(json.dumps(summary, indent=2))
