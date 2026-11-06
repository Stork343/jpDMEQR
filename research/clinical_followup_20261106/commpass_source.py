import json
from pathlib import Path
import subprocess
import sys
from urllib.parse import urlencode

root = Path("data/raw/commpass_20261106")
root.mkdir(parents=True,exist_ok=True)


def gdc(endpoint, filters, fields, size, destination, offset=0):
    url = f"https://api.gdc.cancer.gov/{endpoint}?" + urlencode({"filters":json.dumps(filters),
          "fields":",".join(fields),"size":size,"from":offset,"format":"JSON"})
    if not destination.exists():
        partial = destination.with_suffix(".incoming")
        subprocess.run(["curl","--http1.1","--silent","--show-error","--fail","--location","--retry","5",
                        "--retry-all-errors","--proxy","http://127.0.0.1:7897","--connect-timeout","12","--max-time","100",
                        "--output",str(partial),url],check=True)
        parsed = json.loads(partial.read_text())
        assert "data" in parsed
        partial.replace(destination)
    return json.loads(destination.read_text())


def member(field, values):
    return {"op":"in","content":{"field":field,"value":values}}


if __name__ == "__main__":
    filters = {"op":"and","content":[member("cases.project.project_id",["MMRF-COMMPASS"]),
               member("data_type",["Gene Expression Quantification"]),member("access",["open"])]}
    fields = ["file_id","file_name","file_size","analysis.workflow_type","cases.case_id","cases.submitter_id",
              "cases.samples.sample_id","cases.samples.submitter_id","cases.samples.sample_type","cases.samples.days_to_collection"]
    result = gdc("files",filters,fields,2000,root/"transcript_roster.json")
    assert len(result["data"]["hits"]) == result["data"]["pagination"]["total"]
    clinical = member("project.project_id",["MMRF-COMMPASS"])
    fields = ["case_id","submitter_id","follow_ups.days_to_follow_up","follow_ups.molecular_tests.laboratory_test",
              "follow_ups.molecular_tests.test_units","follow_ups.molecular_tests.days_to_test",
              "follow_ups.molecular_tests.biospecimen_type","follow_ups.molecular_tests.test_analyte_type",
              "diagnoses.treatments.days_to_treatment_start","diagnoses.treatments.therapeutic_agents"]
    schema = gdc("cases",clinical,fields,1,root/"timeline_schema.json")
    print("public RNA files",len(result["data"]["hits"]),"clinical fields",schema.get("warnings",{}),flush=True)
