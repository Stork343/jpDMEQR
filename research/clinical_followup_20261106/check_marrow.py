from pathlib import Path
import json

import duckdb
import numpy as np
import pandas as pd

from panel import study_config

study = study_config("commpass_lightchain")
source = Path(study["source_directory"])
inputs = Path(study["input_directory"])
cohort = pd.read_csv(inputs / "roles.csv").query("role == 'discovery'")
db = duckdb.connect()
db.register("cohort", cohort)
db.read_json(str(source / "discovery_clinical.json"), columns={
    "submitter_id": "VARCHAR", "follow_ups": "STRUCT(follow_up_id VARCHAR, days_to_follow_up INTEGER, "
    "molecular_tests STRUCT(laboratory_test VARCHAR, test_result VARCHAR, test_value DOUBLE, "
    "test_units VARCHAR, biospecimen_type VARCHAR, days_to_test INTEGER)[])[]"}).create_view("registry_cases")
assert db.sql("select count(*) from registry_cases").fetchone()[0] == len(cohort)
assert db.sql("select count(*) from registry_cases anti join cohort on submitter_id=patient").fetchone()[0] == 0
db.sql("""
    create macro period(d) as case when d between 45 and 135 then 1 when d between 136 and 225 then 2
      when d between 226 and 315 then 3 end
""")
reference = db.sql("""
    with tests as (
      select c.submitter_id as patient,v.follow_up_id as visit,v.days_to_follow_up as day,
             coalesce(t.days_to_test,day) as measured,ends_with(t.laboratory_test,', Kappa') as kappa,
             t.test_value * case t.test_units when 'mg/dL' then 10 when 'mg/L' then 1
               when 'mcg/mL' then 1 when 'g/L' then 1000 end as amount
      from registry_cases c,unnest(c.follow_ups) as f(v),unnest(v.molecular_tests) as m(t)
      where t.biospecimen_type='Blood' and t.test_result='Test Value Reported' and t.test_value>=0
        and t.laboratory_test in ('Serum Free Immunoglobulin Light Chain, Kappa','Serum Free Immunoglobulin Light Chain, Lambda')),
    visits as (
      select patient,day,visit,median(amount) filter(where kappa) as k,median(amount) filter(where not kappa) as l
      from tests where day between -30 and 30 and measured<=30 group by all having k is not null and l is not null),
    baseline as (
      select patient,day as base_day,median(k) as k,median(l) as l from visits group by patient,day
      qualify row_number() over(partition by patient order by day desc)=1),
    chain as (
      select *,k/nullif(l,0)>1.65 as is_kappa,case when is_kappa then k else l end as involved from baseline
      where (is_kappa or k/nullif(l,0)<0.26) and involved>=100),
    visits_after as (
      select a.patient,a.day,c.base_day,c.is_kappa,c.involved,period(a.day) as phase,median(a.amount) as future
      from tests a join chain c using(patient) where a.kappa=c.is_kappa and period(a.day)=period(a.measured)
      group by all having future is not null),
    selected as (
      select * from visits_after qualify row_number() over(partition by patient,phase order by abs(day-90*phase),day)=1)
    select *,ln(1+future/involved) as ref_y from selected
""").df()
baseline = pd.read_csv(inputs / "time_audit.csv")
rows = pd.read_csv(inputs / "discovery.csv")
genotype = pd.read_csv(inputs / "features.csv.gz", usecols=["sample_id", "C_proliferation_late", "C_proliferation", "C_involved"])
paired = reference.merge(baseline, left_on=["patient", "phase"], right_on=["patient", "period"],
                         validate="one_to_one", suffixes=("_sql", "_stored")).merge(rows, on="sample_id", validate="one_to_one")
paired = paired.merge(genotype, on="sample_id", validate="one_to_one")
assert len(paired) == len(rows) == len(reference) == len(baseline)
assert paired.day_sql.eq(paired.day_stored).all() and paired.base_day_sql.eq(paired.base_day_stored).all()
assert paired.is_kappa.eq(paired.chain.str.endswith("Kappa")).all()
for stored, calculated in [(paired.y, paired.ref_y), (paired.baseline_mg_l, paired.involved), (paired.future_mg_l, paired.future),
                           (paired.C_involved, np.log1p(paired.involved))]:
    np.testing.assert_allclose(stored, calculated, atol=1e-12, rtol=1e-14)

program = {"MKI67", "TOP2A", "BIRC5", "CCNB1", "UBE2C"}
annotation = pd.read_csv(source / "first_star.tsv", sep="\t", comment="#", usecols=["gene_id", "gene_name", "gene_type"])
names = annotation.loc[annotation.gene_name.isin(program) & annotation.gene_type.eq("protein_coding"), ["gene_id", "gene_name"]]
assert set(names.gene_name) == program
samples = cohort.set_index("patient").loc[paired.patient.unique(), "sample"]
pieces = []
for chunk in pd.read_csv(source / "transcripts.tsv.gz", sep="\t", chunksize=5000, usecols=["Ensembl_ID", *samples]):
    pieces.append(names.merge(chunk, left_on="gene_id", right_on="Ensembl_ID", validate="one_to_one"))
scores = pd.concat(pieces).groupby("gene_name")[list(samples)].median()
assert set(scores.index) == program
baseline_score = paired["sample"].map(scores.mean(skipna=False))
contrast = baseline_score * (paired.phase-1) / 2
np.testing.assert_allclose(paired.C_proliferation, baseline_score, rtol=0, atol=1e-12)
np.testing.assert_allclose(paired.C_proliferation_late, contrast, rtol=0, atol=1e-12)
state = {"source": "independent SQL reconstruction from GDC nested records and a separate five-gene extraction from Xena",
         "patients": paired.patient.nunique(), "records": len(paired),
         "maximum_endpoint_difference": float(abs(paired.y-paired.ref_y).max()),
         "maximum_contrast_difference": float(abs(paired.C_proliferation_late-contrast).max()),
         "baseline_chain_and_units_verified": True, "reserved_outcomes_used": False}
(Path(study["report_directory"]) / "endpoint_verification.json").write_text(json.dumps(state, indent=2) + "\n")
print(json.dumps(state))
