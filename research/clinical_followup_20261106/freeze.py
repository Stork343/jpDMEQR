from pathlib import Path
import csv
import gzip
import json
import subprocess
from zipfile import ZipFile

import numpy as np
import pandas as pd
import yaml

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
    with gzip.open(source / "gene_header.part", "rt") as stream:
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


application = template("application_template")
method = template("method_template")
shared = ["comparators", "target_candidates", "primary_prediction_rule", "working_penalties", "gap_tolerance",
          "elastic_net_alpha", "elastic_net_path_length", "elastic_net_lambda_min_ratio",
          "ridge_variance_ratio_max", "ridge_variance_ratio_min", "forest", "xgboost", "lightgbm"]
for name, definition in registry["questions"].items():
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
               "gaitrec_propulsion": gait_roles, "vitaldb_crp": surgical_roles, "ispy2_volume": breast_roles}
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
