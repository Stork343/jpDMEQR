from pathlib import Path
import pandas as pd
import numpy as np


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
