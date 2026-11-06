from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path

import pandas as pd

from commpass_source import gdc, member
from panel import study_config

study = study_config("commpass_lightchain")
inputs = Path(study["input_directory"])
source = Path(study["source_directory"])
roles = pd.read_csv(inputs/"roles.csv",dtype=str)
cohort = roles.loc[roles.role.eq("discovery")].sort_values("gdc_case")
parts = source/"discovery_clinical_core"
parts.mkdir(exist_ok=True)
fields = ["case_id","submitter_id","demographic.age_at_index","diagnoses.age_at_diagnosis",
          "diagnoses.treatments.therapeutic_agents","diagnoses.treatments.days_to_treatment_start",
          "follow_ups.follow_up_id","follow_ups.days_to_follow_up","follow_ups.molecular_tests.laboratory_test",
          "follow_ups.molecular_tests.test_units","follow_ups.molecular_tests.test_value",
          "follow_ups.molecular_tests.test_result","follow_ups.molecular_tests.days_to_test",
          "follow_ups.molecular_tests.biospecimen_type"]


def download(start):
    ids = cohort.gdc_case.iloc[start:start+12].to_list()
    result = gdc("cases",member("case_id",ids),fields,len(ids),parts/f"batch_{start:04d}.json")
    assert not result.get("warnings"),result.get("warnings")
    records = result["data"]["hits"]
    assert set(r["case_id"] for r in records)==set(ids)
    return records


with ThreadPoolExecutor(max_workers=5) as executor:
    cases = [case for group in executor.map(download,range(0,len(cohort),12)) for case in group]
assert set(c["submitter_id"] for c in cases)==set(cohort.patient)
(source/"discovery_clinical.json").write_text(json.dumps(cases,separators=(",",":")))
ids = {case["submitter_id"] for case in cases}
state = {"public_cases_downloaded":len(ids),"reserved_cases_requested":len(ids & set(roles.loc[roles.role.eq("reserved"),"patient"])),
         "source":"GDC cases endpoint","fields":fields}
(Path(study["report_directory"])/"clinical_source.json").write_text(json.dumps(state,indent=2)+"\n")
print(json.dumps(state))
