import json
from pathlib import Path
import sys

import duckdb
import numpy as np
import pandas as pd

from panel import study_config

name = sys.argv[1]
study = study_config(name)
folder = Path(study["source_directory"]) / "msk_chord_2024"
inputs = Path(study["input_directory"])
roles = pd.read_csv(inputs / "roles.csv").query("role == 'discovery'")
marker = study.get("marker", "psa")
floor = study.get("report_floor", 0.1)
threshold = study.get("baseline_threshold", 2)
unit = "ng/ml" if marker in {"psa", "cea"} else "units/ml"
db = duckdb.connect()
db.register("cohort", roles)
db.read_csv(str(folder / f"data_timeline_{marker}_labs.txt"), delimiter="\t", all_varchar=True).create_view("assays")
wrong_units = db.execute("""
    select count(*) from assays a join cohort c on c.patient=a.PATIENT_ID
    where lower(a.TEST)=? and lower(trim(a.LR_UNIT_MEASURE)) is distinct from ?
""", [marker, unit]).fetchone()[0]
assert wrong_units == 0
measures = db.execute("""
    with typed as (
      select c.patient,c.index_day,c.specimen_day,c.SAMPLE_ID,
             cast(a.START_DATE as integer) as day,trim(a.RESULT) as report
      from assays a join cohort c on c.patient=a.PATIENT_ID where lower(a.TEST)=?),
    amounts as (
      select *,regexp_extract(report,'^(<=|>=|[<≤=≥>])',1) as qualifier,
        try_cast(trim(regexp_replace(report,'^(<=|>=|[<≤=≥>])','')) as double) as amount
      from typed),
    dates as (
      select patient,index_day,specimen_day,SAMPLE_ID,day,day-index_day as elapsed,
             median(greatest(amount,?)) as value from amounts
      where amount is not null and amount>=0 and
        (qualifier in ('','=') or (qualifier in ('<','<=','≤') and amount<=?))
      group by all),
    baseline as (
      select patient,day as base_day,value as baseline from dates where elapsed between -60 and 0
      qualify row_number() over(partition by patient order by day desc)=1),
    windows as (
      select d.*,b.baseline,b.base_day,case when elapsed<=42 then 1 when elapsed<=70 then 2 else 3 end as period
      from dates d join baseline b using(patient) where d.elapsed between 14 and 112 and b.baseline>=?)
    select *,ln(1+value/baseline) as ref_y from windows qualify
      row_number() over(partition by patient,period order by abs(elapsed-28*period),day)=1
""", [marker, floor, floor, threshold]).df()

genes = study.get("marker_genes", ["TP53"])
fields = ["Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification", "HGVSp_Short"]
events = pd.read_csv(folder / "data_mutations.txt", sep="\t", comment="#", dtype=str, usecols=fields)
db.register("events", events.loc[events.Tumor_Sample_Barcode.isin(measures.SAMPLE_ID) & events.Hugo_Symbol.isin(genes)])
protein_altering = db.sql("""
    select * from events where Variant_Classification in ('Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del',
        'In_Frame_Ins','Missense_Mutation','Nonsense_Mutation','Nonstop_Mutation','Splice_Site','Translation_Start_Site')
""")
if marker in {"psa", "ca_19-9"}:
    columns = {"Hugo_Symbol", *measures.SAMPLE_ID}
    cn = pd.read_csv(folder / "data_cna.txt", sep="\t", comment="#",
                     usecols=lambda column: column in columns).set_index("Hugo_Symbol")
    copy = measures.SAMPLE_ID.map(cn.loc[genes[0]])
    altered = measures.SAMPLE_ID.isin(protein_altering.df().Tumor_Sample_Barcode)
    measures["genotype"] = (altered | copy.eq(-2)).astype(float).where(altered | copy.notna())
else:
    sites = protein_altering.filter(
        "HGVSp_Short in ('p.E380Q','p.Y537S','p.Y537N','p.Y537C','p.D538G')" if marker == "ca_15-3" else
        r"regexp_matches(HGVSp_Short,'^p\.[A-Z](12|13|61|117|146)([A-Z*]|del|ins|fs)')").df()
    measures["genotype"] = measures.SAMPLE_ID.isin(sites.Tumor_Sample_Barcode).astype(float)

stored = pd.read_csv(inputs / "discovery.csv")
audit = pd.read_csv(inputs / "time_audit.csv")
x = pd.read_csv(inputs / "features.csv.gz", usecols=["sample_id", study["inference_target"]])
audit = audit.rename(columns={"sample_day": "specimen_day"}) if marker == "psa" else audit
matched = measures.merge(audit, on=["patient", "period"], validate="one_to_one", suffixes=("_sql", "_stored"))
matched = matched.merge(stored, on="sample_id", validate="one_to_one").merge(x, on="sample_id", validate="one_to_one")
assert len(matched) == len(audit) == len(stored) == len(measures)
assert matched.SAMPLE_ID.eq(matched["sample"]).all()
assert matched.day.eq(matched.outcome_day).all() and matched.base_day.eq(matched.baseline_day).all()
assert matched.specimen_day_sql.eq(matched.specimen_day_stored).all() and matched.index_day_sql.eq(matched.index_day_stored).all()
assert matched.at_limit.eq(matched.value.eq(floor)).all()
expected = matched.genotype * (matched.period - 1) / 2
np.testing.assert_allclose(matched.y, matched.ref_y, atol=1e-12, rtol=0)
np.testing.assert_allclose(matched[study["inference_target"]], expected, atol=1e-12, rtol=0, equal_nan=True)
state = {"source": "separate SQL reconstruction of reportable assay records, follow-up windows and genomic contrast",
         "patients": matched.patient.nunique(), "records": len(matched),
         "maximum_endpoint_difference": float((matched.y-matched.ref_y).abs().max()),
         "maximum_contrast_difference": float((matched[study["inference_target"]]-expected).abs().max()),
         "specimen_and_assay_days_verified": True, "reserved_outcomes_used": False}
(Path(study["report_directory"]) / "endpoint_verification.json").write_text(json.dumps(state, indent=2) + "\n")
print(json.dumps(state))
