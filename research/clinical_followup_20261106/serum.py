from pathlib import Path
import json
import pandas as pd
import numpy as np


def genomic_inputs(source, cohort):
    folder = source/"msk_chord_2024"
    covered = []
    for name in ["IMPACT341","IMPACT410","IMPACT468","IMPACT505"]:
        panel = json.loads((source/"gene_panels"/f"{name}.json").read_text())
        assert panel["genePanelId"] == name
        covered.append({g["hugoGeneSymbol"] for g in panel["genes"]})
    genes = sorted(set.intersection(*covered))
    samples = set(cohort.SAMPLE_ID)
    cn = pd.read_csv(folder/"data_cna.txt", sep="\t", comment="#",
                     usecols=lambda x: x=="Hugo_Symbol" or x in samples,index_col="Hugo_Symbol")
    assert cn.index.is_unique
    cn = cn.reindex(genes).T.reindex(cohort.SAMPLE_ID)
    assert np.isin(cn.stack().to_numpy(),[-2,-1.5,-1,0,1,2]).all()
    mutation = pd.DataFrame(0.0,index=cn.index,columns=genes)
    allowed = {"Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del",
               "In_Frame_Ins","Splice_Site","Nonstop_Mutation","Translation_Start_Site"}
    records = []
    for chunk in pd.read_csv(folder/"data_mutations.txt",sep="\t",comment="#",dtype=str,chunksize=50000,
                             usecols=["Hugo_Symbol","Tumor_Sample_Barcode","Variant_Classification","HGVSp_Short"]):
        q = chunk.loc[chunk.Tumor_Sample_Barcode.isin(samples) & chunk.Hugo_Symbol.isin(genes)
                      & chunk.Variant_Classification.isin(allowed)]
        records.append(q)
    events = pd.concat(records,ignore_index=True)
    for sample,gene in events[["Tumor_Sample_Barcode","Hugo_Symbol"]].drop_duplicates().itertuples(index=False,name=None):
        mutation.at[sample,gene] = 1
    cn.index = cohort.index
    mutation.index = cohort.index
    events["patient"] = events.Tumor_Sample_Barcode.map(dict(zip(cohort.SAMPLE_ID,cohort.index)))
    return cn,mutation,events


def marker_series(folder, cohort, marker, floor, unit):
    cols = ["PATIENT_ID","START_DATE","TEST","RESULT","LR_UNIT_MEASURE"]
    raw = pd.read_csv(folder/f"data_timeline_{marker}_labs.txt",sep="\t",comment="#",dtype=str,usecols=cols)
    selected = raw.loc[raw.PATIENT_ID.isin(cohort.index) & raw.TEST.str.lower().eq(marker)].copy()
    assert selected.LR_UNIT_MEASURE.str.strip().str.lower().eq(unit).all()
    selected["day"] = pd.to_numeric(selected.START_DATE)
    groups = selected.RESULT.str.extract(r"^\s*(<=|<|≤|=|>=|>|≥)?\s*(\d+(?:\.\d*)?|\.\d+)\s*$")
    numbers = pd.to_numeric(groups[1])
    exact = groups[0].isna() | groups[0].eq("=")
    bounded = groups[0].isin(["<","<=","≤"]) & numbers.le(floor)
    selected["value"] = numbers.clip(lower=floor).where(exact|bounded)
    count = int(selected.value.isna().sum())
    frame = selected.dropna(subset=["day","value"]).groupby(["PATIENT_ID","day"],as_index=False).value.median()
    frame["elapsed"] = frame.day-frame.PATIENT_ID.map(cohort.index_day)
    return frame,count


def serum_roles(source, seed, study, reserved_fraction, number_folds):
    folder = source / "msk_chord_2024"
    def read(table, columns):
        return pd.read_csv(folder/("data_"+table+".txt"),sep="\t",comment="#",dtype=str,usecols=columns)
    specimens = read("clinical_sample", ["PATIENT_ID","SAMPLE_ID","CANCER_TYPE","CANCER_TYPE_DETAILED","SAMPLE_TYPE","GENE_PANEL"])
    dates = read("timeline_specimen", ["PATIENT_ID","SAMPLE_ID","START_DATE"])
    specimens = specimens.merge(dates,on=["PATIENT_ID","SAMPLE_ID"], validate="one_to_one")
    specimens["specimen_day"] = pd.to_numeric(specimens.pop("START_DATE"))
    first = specimens.loc[specimens.specimen_day.eq(specimens.groupby("PATIENT_ID").specimen_day.transform("min"))]
    first = first.loc[first.groupby("PATIENT_ID").CANCER_TYPE.transform("nunique").eq(1),["PATIENT_ID","CANCER_TYPE"]].drop_duplicates()
    specimens = specimens.merge(first,on=["PATIENT_ID","CANCER_TYPE"],validate="many_to_one")
    prostate = pd.read_csv("data/derived/chord_psa_20261106/roles.csv",usecols=["patient"],dtype=str)
    specimens = specimens.loc[~specimens.PATIENT_ID.isin(prostate.patient)
        & specimens.CANCER_TYPE_DETAILED.isin(study["histologies"])
        & specimens.GENE_PANEL.isin(["IMPACT341","IMPACT410","IMPACT468","IMPACT505"])]
    profiles = read("gene_panel_matrix",["SAMPLE_ID","mutations","cna"])
    specimens = specimens.merge(profiles,on="SAMPLE_ID",validate="one_to_one")
    specimens = specimens.loc[specimens.mutations.eq(specimens.GENE_PANEL) & specimens.cna.eq(specimens.GENE_PANEL)]
    therapy = read("timeline_treatment",["PATIENT_ID","START_DATE","AGENT"])
    therapy = therapy.loc[therapy.AGENT.isin(study["index_agents"]) & therapy.START_DATE.notna()].copy()
    therapy["index_day"] = pd.to_numeric(therapy.START_DATE)
    starts = therapy.groupby(["PATIENT_ID","index_day"]).AGENT.agg(lambda s:"+".join(sorted(set(s)))).reset_index()
    gap = starts.groupby("PATIENT_ID").index_day.diff()
    starts = starts.loc[gap.isna()|gap.gt(90)]
    frame = starts.merge(specimens,on="PATIENT_ID")
    frame = frame.loc[(frame.index_day-frame.specimen_day).between(30,730)].sort_values(
        ["PATIENT_ID","index_day","specimen_day","SAMPLE_ID"],ascending=[True,True,False,True]).drop_duplicates("PATIENT_ID")
    frame = frame.rename(columns={"PATIENT_ID":"patient"})
    frame["role"] = "discovery"
    frame["outer_fold"] = 0
    rng = np.random.default_rng(seed)
    for _, cell in frame.groupby(["AGENT",frame.SAMPLE_TYPE.eq("Metastasis")]):
        ids = rng.permutation(cell.sort_values("patient").index)
        n = round(len(ids)*reserved_fraction)
        frame.loc[ids[:n],"role"] = "reserved"
        frame.loc[ids[n:],"outer_fold"] = (np.arange(len(ids)-n)+rng.integers(number_folds))%number_folds+1
    return frame.sort_values("patient")
