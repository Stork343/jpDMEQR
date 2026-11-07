from concurrent.futures import ThreadPoolExecutor
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import yaml

from panel import settings_directory, study_config

parser = argparse.ArgumentParser()
parser.add_argument("study")
parser.add_argument("--implementation",type=Path)
parser.add_argument("--workers",type=int,default=2,choices=[1,2,3,4])
args = parser.parse_args()
name = args.study
study = study_config(name)
root = Path(study["output_root"])
method = yaml.safe_load((settings_directory/f"cq_{name}.yml").read_text())
implementation = args.implementation or root/"implementation"
assert yaml.safe_load((implementation/"config.yml").read_text()) == study
jobs = pd.read_csv(root/"jobs.csv").query("method == 'JOINT_DIRECT'")
jobs = jobs.assign(priority=jobs.tau.ne(study["decision_quantile"])).sort_values(["priority","trait","tau"])
assert len(jobs) == study["outer_folds"]*len(study["quantiles"])
env = dict(os.environ,TMPDIR=str(Path("results/scratch/runtime").resolve()),OMP_NUM_THREADS="1",
           OPENBLAS_NUM_THREADS="1",VECLIB_MAXIMUM_THREADS="1",R_MAX_VSIZE="24Gb")


def fit(job):
    directory = Path(job.directory)
    if not (directory/"status.json").exists():
        assert not directory.exists(), str(directory)
        with (root/f"{job.trait}_JOINT_DIRECT_{job.run_id}.log").open("a") as stream:
            subprocess.run(["Rscript",str(implementation/"run_discovery.R"),job.trait,job.method,
                            job.run_id,str(job.tau),str(implementation.resolve())],check=True,env=env,
                           stdout=stream,stderr=subprocess.STDOUT)
    assert json.loads((directory/"status.json").read_text())["status"] == "completed"
    print(job.trait,job.run_id,"completed",flush=True)
    return pd.read_csv(directory/"prediction.csv")


with ThreadPoolExecutor(max_workers=args.workers) as workers:
    predictions = pd.concat(workers.map(fit,jobs.itertuples()),ignore_index=True)
new = pd.read_csv(Path(method["output_root"])/study["study_id"]/"predictions.csv")
rows = pd.read_csv(Path(study["input_directory"])/"discovery.csv")
assert rows.role.eq("discovery").all()
records = pd.concat([predictions,new[["sample_id","strain_id","method","tau","prediction"]]],ignore_index=True)
assert not records.duplicated(["sample_id","method","tau"]).any()
wide = records.pivot(index=["method","sample_id","strain_id"],columns="tau",values="prediction")
assert list(wide.columns) == study["quantiles"] and wide.notna().all().all()
q = np.maximum(0,np.sort(wide.to_numpy(),axis=1))
value = pd.DataFrame({"prediction":q[:,study["quantiles"].index(study["decision_quantile"])]},index=wide.index).reset_index()
value = value.merge(rows[["sample_id","strain_id","y"]],on=["sample_id","strain_id"],validate="many_to_one")
residual = value.y-value.prediction
tau = study["decision_quantile"]
value["loss"] = np.maximum(tau*residual,(tau-1)*residual)
patient = value.groupby(["method","strain_id"]).loss.mean().unstack("method")
assert not patient.isna().any().any() and len(patient) == rows.strain_id.nunique()
loss = patient.mean()
gain = 1-loss[method["prediction_method"]]/loss.JOINT_DIRECT
state = {"study_id":study["study_id"],"patients":len(patient),"records":len(rows),"quantile":tau,
         "cq_loss":float(loss[method["prediction_method"]]),"joint_direct_loss":float(loss.JOINT_DIRECT),
         "relative_gain":float(gain),"passed":bool(gain>=study["screening_minimum_relative_gain"]),
         "formal_validation_used":False,"manuscript_adoption":False}
report = Path(study["report_directory"])
progress = jobs[["trait","method","directory"]].copy()
progress["status"] = "completed"
progress["exit_code"] = 0
progress.to_csv(root/"joint_direct_progress.csv",index=False)
loss_table = []
for scoring, frame in [("raw", records),("primary",pd.DataFrame(q,index=wide.index,columns=wide.columns).stack().rename("prediction").reset_index())]:
    frame = frame.merge(rows[["sample_id","strain_id","y"]],on=["sample_id","strain_id"],validate="many_to_one")
    error = frame.y-frame.prediction
    frame["loss"] = np.maximum(frame.tau*error,(frame.tau-1)*error)
    group = frame.groupby(["method","tau","strain_id"]).loss.mean().reset_index()
    stats = group.groupby(["method","tau"]).loss.agg(patients="size",loss="mean",sd="std").reset_index()
    stats["se"] = stats.sd/np.sqrt(stats.patients)
    stats["scoring"] = scoring
    loss_table.append(stats)
pd.concat(loss_table,ignore_index=True).to_csv(report/"joint_direct_loss.csv",index=False)
subprocess.run(["Rscript","research/clinical_search_20261105/verify_screen.R",
                str(settings_directory/f"{name}.yml"),"joint_direct_progress.csv","joint_direct_loss.csv"],
               check=True,env=env)
assert json.loads((report/"joint_direct_screen_audit.json").read_text())["status"] == "passed"
(report/"joint_direct_comparison.json").write_text(json.dumps(state,indent=2)+"\n")
print(json.dumps(state))
