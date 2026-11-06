from pathlib import Path
import csv
import gzip
import json
import subprocess
import sys
from zipfile import ZipFile

import numpy as np
import pandas as pd
import yaml
from serum import serum_roles

registry = yaml.safe_load(Path("config/clinical_followup_search_20261106.yml").read_text())
root = Path("results/raw/clinical_followup_settings_20261106")
root.mkdir(parents=True, exist_ok=True)


def template(name):
    specification = f"{registry['template_commit']}:{registry[name]}"
    return yaml.safe_load(subprocess.check_output(["git", "show", specification]))


def immutable(path, content):
    if path.exists():
        assert path.read_text() == content, str(path)
    else:
        path.write_text(content)


def allocate_sites(counts, weights, random):
    order = list(random.permutation(sorted(counts.index)))
    size = np.zeros(len(weights))
    assignment = {}
    for site in sorted(order, key=lambda s: -counts[s]):
        group = int(np.argmin(size / np.asarray(weights)))
        assignment[site] = group
        size[group] += counts[site]
    return assignment


def eye_roles(source, seed):
    with ZipFile(source / "release.zip") as archive:
        member, = [n for n in archive.namelist() if n.endswith("/CSV/VF_Data.csv")]
        with archive.open(member) as stream:
            ids = pd.read_csv(stream, usecols=["PatID"], dtype=str).PatID.unique()
    random = np.random.default_rng(seed)
    shuffled = random.permutation(sorted(ids))
    hold = round(len(shuffled) * registry["validation_fraction"])
    frame = pd.DataFrame({"patient": shuffled, "role": "discovery", "outer_fold": 0})
    frame.loc[:hold - 1, "role"] = "reserved"
    frame.loc[hold:, "outer_fold"] = np.arange(len(frame) - hold) % registry["outer_folds"] + 1
    return frame.sort_values("patient")


def icu_roles(source, seed):
    cols = ["uniquepid", "patientunitstayid", "patienthealthsystemstayid", "unitvisitnumber",
            "hospitaldischargeyear", "hospitalid", "age", "gender"]
    roster = pd.read_csv(source / "patient.csv.gz", usecols=cols, dtype={"uniquepid": str})
    roster["age"] = pd.to_numeric(roster.age.replace({"> 89": 90, ">89": 90}), errors="raise")
    roster = roster.loc[roster.age.ge(18)].sort_values(
        ["hospitaldischargeyear", "patienthealthsystemstayid", "unitvisitnumber", "patientunitstayid"])
    frame = roster.drop_duplicates("uniquepid").rename(columns={"uniquepid": "patient"}).copy()
    assert frame.patient.notna().all() and frame.patientunitstayid.is_unique
    random = np.random.default_rng(seed)
    counts = frame.hospitalid.value_counts()
    reserved = allocate_sites(counts, [1 - registry["validation_fraction"], registry["validation_fraction"]], random)
    frame["role"] = frame.hospitalid.map(reserved).map({0: "discovery", 1: "reserved"})
    folds = allocate_sites(frame.loc[frame.role.eq("discovery"), "hospitalid"].value_counts(),
                           [1] * registry["outer_folds"], random)
    frame["outer_fold"] = frame.hospitalid.map(folds).fillna(-1).astype(int) + 1
    assert frame.groupby("hospitalid").role.nunique().eq(1).all()
    assert frame.groupby("hospitalid").outer_fold.nunique().eq(1).all()
    return frame.sort_values("patient")


def gait_roles(source, seed):
    metadata = pd.read_csv(source / "GRF_metadata.csv", dtype={"SUBJECT_ID": str})
    patients = metadata.loc[metadata.CLASS_LABEL.isin(["A", "C", "H", "K"])].copy()
    assert patients.TRAIN.add(patients.TEST).eq(1).all()
    assert patients.groupby("SUBJECT_ID").TEST.nunique().eq(1).all()
    patients["date"] = pd.to_datetime(patients.SESSION_DATE, format="%d.%m.%Y")
    frame = patients.sort_values(["date", "SESSION_ID"]).drop_duplicates("SUBJECT_ID")
    frame = frame[["SUBJECT_ID", "CLASS_LABEL", "TEST"]].rename(columns={"SUBJECT_ID": "patient"}).copy()
    frame["role"] = frame.TEST.map({0: "discovery", 1: "reserved"})
    frame["outer_fold"] = 0
    rng = np.random.default_rng(seed)
    for _, group in frame.loc[frame.role.eq("discovery")].groupby("CLASS_LABEL"):
        indices = rng.permutation(group.sort_values("patient").index)
        frame.loc[indices, "outer_fold"] = np.arange(len(indices)) % registry["outer_folds"] + 1
    return frame.drop(columns="TEST").sort_values("patient")


def surgical_roles(source, seed):
    columns = ["caseid", "subjectid", "department", "optype", "approach", "age", "opstart", "opend", "caseend", "aneend"]
    cases = pd.read_csv(source / "cases.csv.gz", usecols=columns, dtype={"subjectid": str})
    types = ["Colorectal", "Biliary/Pancreas", "Stomach", "Hepatic"]
    cases = cases.loc[cases.department.eq("General surgery") & cases.age.ge(18) & cases.optype.isin(types)
                      & cases.approach.isin(["Open", "Videoscopic", "Robotic"])
                      & cases.opend.gt(cases.opstart) & cases.caseend.ge(cases.opend)].copy()
    cases["landmark"] = cases[["caseend", "aneend", "opend"]].max(axis=1)
    cases = cases.sort_values("caseid").drop_duplicates("subjectid").rename(columns={"subjectid": "patient"})
    cases["role"] = "discovery"
    cases["outer_fold"] = 0
    random = np.random.default_rng(seed)
    for _, cell in cases.groupby(["optype", "approach"]):
        index = random.permutation(cell.sort_values("patient").index)
        n = round(len(cell) * registry["validation_fraction"])
        cases.loc[index[:n], "role"] = "reserved"
        cases.loc[index[n:], "outer_fold"] = (np.arange(len(index)-n) + random.integers(5)) % 5 + 1
    assert cases.caseid.is_unique and cases.patient.notna().all()
    return cases.sort_values("patient")


def breast_roles(source, seed):
    clinical = pd.read_excel(source / "ISPY2-Imaging-Cohort-1-Clinical-Data.xlsx", dtype={"Patient_ID": str},
                             usecols=["Patient_ID", "HR", "HER2", "Arm", "MP", "Age_at_Screening"])
    ids = pd.read_excel(source / "Multi-feature-MRI-NACT-Data.xlsx", dtype=str, usecols=["CLINICAL-TRIAL-SUBJECT-ID"])
    gene_source = source / "expression.tsv.gz" if (source / "expression.tsv.gz").exists() else source / "gene_header.part"
    with gzip.open(gene_source, "rt") as stream:
        genes = next(csv.reader([stream.readline()], delimiter="\t"))
    assert clinical.Patient_ID.is_unique and ids.iloc[:, 0].is_unique and len(set(genes)) == len(genes)
    frame = clinical.loc[clinical.Patient_ID.isin(set(genes) & set(ids.iloc[:, 0]))].rename(columns={"Patient_ID": "patient"}).copy()
    frame["role"] = "discovery"
    frame["outer_fold"] = 0
    random = np.random.default_rng(seed)
    for _, cell in frame.groupby(["HR", "HER2", "Arm"]):
        index = random.permutation(cell.sort_values("patient").index)
        n = round(len(index) * registry["validation_fraction"])
        frame.loc[index[:n], "role"] = "reserved"
        frame.loc[index[n:], "outer_fold"] = (np.arange(len(index)-n) + random.integers(5)) % 5 + 1
    return frame.sort_values("patient")


def protein_roles(source, seed):
    original = breast_roles(source, seed)
    protein = pd.read_csv(source / "signaling_roster.csv", dtype={"patient": str})
    assert protein.patient.is_unique and len(protein) == 736
    return original.merge(protein, on="patient", how="inner", validate="one_to_one").sort_values("patient")


def registry_roles(source, seed):
    cols = ["PtID","SiteID","Visit","ExamDtMnC","ageConsent","diabDur","BGTestAvgNumMeter","BGTestAvgNumPtRep"]
    with ZipFile(source / "registry.zip") as archive:
        name, = [s for s in archive.namelist() if s.endswith("/VisitsLong.txt")]
        frame = pd.read_csv(archive.open(name), sep="|", usecols=cols, dtype={"PtID": str})
    frame = frame.loc[frame.Visit.isin(["Enrollment1","Enrollment2"]) & frame.ageConsent.between(14,24)
                      & frame.ExamDtMnC.between(-1,1) & frame.diabDur.ge(1)
                      & frame.BGTestAvgNumMeter.between(0,30) & frame.BGTestAvgNumPtRep.between(0,30)
                      & frame.SiteID.notna()].copy()
    frame = frame.sort_values(["PtID","ExamDtMnC","Visit"]).drop_duplicates("PtID").rename(columns={"PtID":"patient"})
    rng = np.random.default_rng(seed)
    sites = allocate_sites(frame.SiteID.value_counts(), [1-registry["validation_fraction"], registry["validation_fraction"]], rng)
    frame["role"] = frame.SiteID.map(sites).map({0:"discovery",1:"reserved"})
    groups = allocate_sites(frame.loc[frame.role.eq("discovery"),"SiteID"].value_counts(), [1]*registry["outer_folds"], rng)
    frame["outer_fold"] = frame.SiteID.map(groups).fillna(-1).astype(int)+1
    assert frame.groupby("SiteID").role.nunique().eq(1).all() and frame.groupby("SiteID").outer_fold.nunique().eq(1).all()
    return frame.sort_values("patient")


def prostate_roles(source, seed):
    folder = source / "msk_chord_2024"
    def read(file, columns):
        return pd.read_csv(folder / ("data_"+file+".txt"), sep="\t", comment="#", usecols=columns, dtype=str)
    sample = read("clinical_sample", ["PATIENT_ID","SAMPLE_ID","SAMPLE_TYPE","GENE_PANEL","CANCER_TYPE_DETAILED"])
    sample = sample.loc[sample.CANCER_TYPE_DETAILED.eq("Prostate Adenocarcinoma")
                        & sample.GENE_PANEL.isin(["IMPACT341","IMPACT410","IMPACT468","IMPACT505"])]
    panel = read("gene_panel_matrix", ["SAMPLE_ID","mutations","cna"])
    sample = sample.merge(panel, on="SAMPLE_ID", validate="one_to_one")
    sample = sample.loc[sample.mutations.eq(sample.GENE_PANEL) & sample.cna.eq(sample.GENE_PANEL)]
    time = read("timeline_specimen", ["PATIENT_ID","SAMPLE_ID","START_DATE"])
    sample = sample.merge(time, on=["PATIENT_ID","SAMPLE_ID"],validate="one_to_one")
    sample["specimen_day"] = pd.to_numeric(sample.pop("START_DATE"))
    treatment = read("timeline_treatment", ["PATIENT_ID","START_DATE","AGENT"])
    treatment["index_day"] = pd.to_numeric(treatment.pop("START_DATE"))
    agents = ["ABIRATERONE","ENZALUTAMIDE","APALUTAMIDE","DAROLUTAMIDE"]
    treatment = treatment.loc[treatment.AGENT.isin(agents) & treatment.index_day.notna()].drop_duplicates()
    first = treatment.loc[treatment.index_day.eq(treatment.groupby("PATIENT_ID").index_day.transform("min"))].copy()
    first = first.loc[first.groupby("PATIENT_ID").AGENT.transform("nunique").eq(1)].drop_duplicates("PATIENT_ID")
    frame = first.merge(sample, on="PATIENT_ID")
    frame = frame.loc[(frame.index_day-frame.specimen_day).between(30,730)].sort_values(
        ["PATIENT_ID","specimen_day","SAMPLE_ID"],ascending=[True,False,True]).drop_duplicates("PATIENT_ID")
    frame = frame.rename(columns={"PATIENT_ID":"patient"})
    frame["role"] = "discovery"
    frame["outer_fold"] = 0
    rng = np.random.default_rng(seed)
    for _, cell in frame.groupby(["AGENT",frame.SAMPLE_TYPE.eq("Metastasis")]):
        index = rng.permutation(cell.sort_values("patient").index)
        n = round(len(index)*registry["validation_fraction"])
        frame.loc[index[:n],"role"] = "reserved"
        frame.loc[index[n:],"outer_fold"] = (np.arange(len(index)-n)+rng.integers(5))%5+1
    return frame.sort_values("patient")


def myeloma_roles(source, seed):
    records = json.loads((source/"transcript_roster.json").read_text())["data"]["hits"]
    with gzip.open(source/"transcripts.tsv.gz","rt") as stream:
        samples = next(csv.reader([stream.readline()],delimiter="\t"))[1:]
    assert len(set(samples)) == len(samples)
    rows = []
    for record in records:
        for case in record["cases"]:
            for sample in case.get("samples",[]):
                if (record["analysis"]["workflow_type"] == "STAR - Counts" and
                    sample["sample_type"] == "Primary Blood Derived Cancer - Bone Marrow" and
                    sample["submitter_id"].endswith("_1_BM_CD138pos") and sample["submitter_id"] in samples):
                    rows.append((case["submitter_id"],case["case_id"],sample["submitter_id"]))
    frame = pd.DataFrame(rows,columns=["patient","gdc_case","sample"]).drop_duplicates().sort_values("patient").reset_index(drop=True)
    assert frame.patient.is_unique and frame.gdc_case.is_unique and frame["sample"].is_unique
    rng = np.random.default_rng(seed)
    order = rng.permutation(frame.index)
    n = round(len(frame)*registry["validation_fraction"])
    frame["role"] = "discovery"
    frame["outer_fold"] = 0
    frame.loc[order[:n],"role"] = "reserved"
    frame.loc[order[n:],"outer_fold"] = np.arange(len(order)-n)%registry["outer_folds"]+1
    return frame


def diameter_roles(source, seed):
    fields = ["SUBJECTID", "HR Pos", "Her2MostPos"]
    trial = pd.read_excel(source / "ispy1_clinical.xlsx", sheet_name="TCIA Patient Clinical Subset",
                         usecols=fields, dtype={"SUBJECTID": str}).rename(columns={"SUBJECTID": "patient"})
    expression = pd.read_csv(source / "ispy1_expression_roster.csv", dtype=str)
    training = expression.merge(trial, on="patient", validate="one_to_one")
    training["source_id"] = training.patient
    training["patient"] = "IS1_" + training.patient
    training["role"] = "discovery"
    training["outer_fold"] = 0
    generator = np.random.default_rng(seed)
    strata = training.groupby(["HR Pos", "Her2MostPos", "platform"], dropna=False)
    for _, cell in strata:
        order = generator.permutation(cell.sort_values("patient").index)
        training.loc[order, "outer_fold"] = (np.arange(len(order))+generator.integers(5)) % 5 + 1
    held = pd.read_csv("data/derived/ispy2_volume_20261106/roles.csv", dtype={"patient": str})
    held = held.loc[held.role.eq("reserved"), ["patient", "role", "outer_fold"]].copy()
    assert len(held) == 132 and held.outer_fold.eq(0).all()
    held["source_id"] = held.patient
    held["patient"] = "IS2_" + held.patient
    columns = ["patient", "source_id", "role", "outer_fold"]
    return pd.concat([training[columns], held[columns]]).sort_values("patient")


def acid_base_roles(source, seed):
    catalog = pd.read_csv(source / "catalog.csv", usecols=["patient", "site"])
    assert catalog.patient.is_unique and catalog.groupby("site").size().to_dict() == {"A": 20336, "B": 20000}
    catalog["role"] = catalog.site.map({"A": "discovery", "B": "reserved"})
    catalog["outer_fold"] = 0
    rng = np.random.default_rng(seed)
    order = rng.permutation(catalog.loc[catalog.site.eq("A")].sort_values("patient").index)
    catalog.loc[order, "outer_fold"] = np.arange(len(order)) % registry["outer_folds"] + 1
    return catalog.sort_values("patient")


application = template("application_template")
method = template("method_template")
shared = ["comparators", "target_candidates", "primary_prediction_rule", "working_penalties", "gap_tolerance",
          "elastic_net_alpha", "elastic_net_path_length", "elastic_net_lambda_min_ratio",
          "ridge_variance_ratio_max", "ridge_variance_ratio_min", "forest", "xgboost", "lightgbm"]
requested = sys.argv[1:] or list(registry["questions"])
assert set(requested) <= set(registry["questions"])
for name in requested:
    definition = registry["questions"][name]
    study = {key: application[key] for key in shared}
    study.update(definition)
    study_id = study["study_id"]
    inputs = Path("data/derived") / study_id
    report = Path("results/preflight/empirical_selection_20261101") / study_id
    for directory in [inputs, report]:
        directory.mkdir(parents=True, exist_ok=True)
    study.update(input_directory=str(inputs), output_root=f"results/raw/application/{study_id}",
                 report_directory=str(report), outer_folds=registry["outer_folds"], inner_folds=registry["inner_folds"],
                 quantiles=registry["quantiles"], minimum_observed_fraction=definition.get("minimum_observed_fraction", registry["eligible_fraction"]),
                 screening_minimum_relative_gain=registry["required_relative_gain"], formal_validation_used=False,
                 manuscript_adoption=False)
    study_file = root / f"{name}.yml"
    immutable(study_file, yaml.safe_dump(study, sort_keys=False))
    fit = dict(method, method_id=f"cq_{study_id}", output_root=f"results/raw/cq_{study_id}",
               report_directory=f"results/preflight/cq_{study_id}", source_studies=[str(study_file)],
               inference_target=study["inference_target"], inference_quantile=study["decision_quantile"])
    immutable(root / f"cq_{name}.yml", yaml.safe_dump(fit, sort_keys=False))
    source = Path(study["source_directory"])
    readers = {"uwhvf_central": eye_roles, "eicu_lactate": icu_roles, "eicu_creatinine": icu_roles,
               "gaitrec_propulsion": gait_roles, "vitaldb_crp": surgical_roles, "ispy2_volume": breast_roles,
               "ispy2_signaling": protein_roles, "t1dx_monitoring": registry_roles, "chord_psa": prostate_roles,
               "commpass_lightchain": myeloma_roles, "ispy_diameter": diameter_roles,
               "physionet_ph": acid_base_roles}
    if "marker" in study:
        frame = serum_roles(source, study["role_seed"], study, registry["validation_fraction"], registry["outer_folds"])
    else:
        frame = readers[name](source, study["role_seed"])
    assert frame.patient.is_unique and frame.outer_fold.eq(0).equals(frame.role.eq("reserved"))
    immutable(inputs / "roles.csv", frame.to_csv(index=False))
    inventory = frame.groupby(["role", "outer_fold"]).size().rename("patients")
    immutable(report / "role_inventory.csv", inventory.to_csv())
    state = {"study_id": study_id, "roster_patients": len(frame), "reserved_patients": int(frame.role.eq("reserved").sum()),
             "discovery_patients": int(frame.role.eq("discovery").sum()), "reserved_outcomes_used": False}
    if "hospitalid" in frame:
        state["hospital_codes"] = int(frame.hospitalid.nunique())
        state["reserved_hospital_codes"] = int(frame.loc[frame.role.eq("reserved"), "hospitalid"].nunique())
    immutable(report / "roles.json", json.dumps(state, indent=2) + "\n")
    print(json.dumps(state))
