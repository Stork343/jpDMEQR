import json
from pathlib import Path
import re

import numpy as np
import pandas as pd

from panel import finish_panel, study_config

study = study_config("ispy2_signaling")
inputs = Path(study["input_directory"])
original = study_config("ispy2_volume")
original_inputs = Path(original["input_directory"])
roles = pd.read_csv(inputs / "roles.csv", dtype={"patient": str}).set_index("patient")
reference_roles = pd.read_csv(original_inputs / "roles.csv", dtype={"patient": str}).set_index("patient")
assert roles[["role", "outer_fold"]].equals(reference_roles.loc[roles.index, ["role", "outer_fold"]])
signals = pd.read_csv(inputs / "discovery_signals.csv", dtype={"patient": str}).set_index("patient")
assert signals.index.is_unique and set(signals.index) == set(roles.index[roles.role.eq("discovery")])
assert signals.shape[1] == 139
arm = roles.loc[signals.index, "Arm"]
members = {"AKT": arm.str.contains("MK-2206", regex=False, na=False),
           "HER": arm.str.contains("Neratinib|Pertuzumab|Trastuzumab|T-DM1", regex=True, na=False),
           "PDL1": arm.str.contains("Pembrolizumab", regex=False, na=False)}
target_patterns = {"AKT": r"^(AKT|mTOR|PRAS40|S6RP|p70S6K|X4EBP1)\.[ST][0-9]",
                   "HER": r"^(EGFR|ERBB2|ERBB3)\.Y[0-9]", "PDL1": r"^PDL1\."}
features = {}
catalog = []
for i, antibody in enumerate(signals.columns):
    name = f"P_RPPA_{i:03d}"
    features[name] = signals[antibody]
    catalog.append({"feature": name, "antibody": antibody, "therapy": "all"})
    for group, pattern in target_patterns.items():
        if re.match(pattern, antibody):
            interaction = f"P_{group}_{i:03d}"
            features[interaction] = signals[antibody].where(members[group], 0.0)
            catalog.append({"feature": interaction, "antibody": antibody, "therapy": group})
catalog = pd.DataFrame(catalog)
assert set(catalog.therapy) == {"all", *members}
catalog.to_csv(inputs / "assay_catalog.csv", index=False)
abundance = pd.DataFrame(features)
rows = pd.read_csv(original_inputs / "discovery.csv")
audit = pd.read_csv(original_inputs / "time_audit.csv", dtype={"patient": str})
audit = audit.loc[audit.patient.isin(signals.index)].copy()
rows = rows.merge(audit[["sample_id", "patient"]], on="sample_id", validate="one_to_one")
clinical = pd.read_csv(original_inputs / "features.csv.gz", usecols=lambda col: col == "sample_id" or col.startswith("C_"))
x = rows[["patient", "sample_id"]].merge(clinical, on="sample_id", validate="one_to_one").join(abundance, on="patient")
verification = json.loads((Path(original["report_directory"]) / "endpoint_verification.json").read_text())
assert verification["max_endpoint_difference"] < 1e-10 and verification["max_clinical_input_difference"] < 1e-10
np.testing.assert_allclose(rows.y, np.log1p(audit.set_index("sample_id").loc[rows.sample_id,"followup_volume"].to_numpy() /
    audit.set_index("sample_id").loc[rows.sample_id,"baseline_volume"].to_numpy()), atol=1e-12)
finish_panel(study, rows[["sample_id","patient","y"]], x.drop(columns="patient"), audit)
state = {"source": "published GEO ExpressionSet and verified MRI and gene-expression source match",
         "patients": rows.patient.nunique(), "records": len(rows), "protein_measurements": signals.shape[1],
         "treatment_interactions": len(features)-signals.shape[1], "patient_roles_preserved": True,
         "reserved_outcomes_used": False}
(Path(study["report_directory"]) / "endpoint_verification.json").write_text(json.dumps(state, indent=2)+"\n")
print(json.dumps(state))
