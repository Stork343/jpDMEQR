import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linprog
from scipy.sparse import bmat, csr_matrix, eye
import yaml

from panel import registry, settings_directory

records = []
for name, definition in registry["questions"].items():
    study_id = definition["study_id"]
    weighted = settings_directory/f"pw_{study_id}.yml"
    key = ("labels_t1dx_monitoring" if name == "t1dx_monitoring" else
           ("pw_"+study_id if weighted.exists() else name))
    config_file = settings_directory / f"{key}.yml"
    if not config_file.exists():
        continue
    cfg = yaml.safe_load(config_file.read_text())
    method = yaml.safe_load((settings_directory/f"cq_{key}.yml").read_text())
    report = Path(cfg["report_directory"])
    if not (report/"decision.json").exists():
        continue
    rows = pd.read_csv(Path(cfg["input_directory"])/"discovery.csv").set_index("sample_id")
    assert rows.role.eq("discovery").all()
    tau = cfg["decision_quantile"]
    experts = {}
    for learner in ["CLINICAL_QR","ENET","CQ_SMOOTH"]:
        if learner == "CQ_SMOOTH":
            long = pd.read_csv(Path(method["output_root"])/study_id/"predictions.csv")
        else:
            root = Path(cfg["output_root"])
            long = pd.concat([pd.read_csv(root/f"outer_{i}"/"models"/(learner+"_all")/"prediction.csv")
                              for i in range(1,cfg["outer_folds"]+1)])
        wide = long.pivot(index="sample_id",columns="tau",values="prediction").reindex(rows.index)
        assert np.isfinite(wide).all().all() and list(wide.columns)==cfg["quantiles"]
        experts[learner] = np.maximum(0,np.sort(wide.to_numpy(),axis=1))[:,list(wide.columns).index(tau)]
    matrix = pd.DataFrame(experts,index=rows.index).to_numpy()
    response = rows.y.to_numpy()
    patient_weight = (1/rows.groupby("strain_id").strain_id.transform("size")/rows.strain_id.nunique()).to_numpy()
    n,k = matrix.shape
    constraints = bmat([[csr_matrix(matrix),eye(n),-eye(n)],[csr_matrix(np.ones((1,k))),csr_matrix((1,n)),csr_matrix((1,n))]],
                       format="csr")
    fit = linprog(np.r_[np.zeros(k),tau*patient_weight,(1-tau)*patient_weight], A_eq=constraints,
                  b_eq=np.r_[response,1],bounds=(0,None),method="highs")
    assert fit.success and abs(fit.x[:k].sum()-1)<1e-7
    error = response[:,None] - matrix
    risks = np.sum(np.maximum(tau*error,(tau-1)*error)*patient_weight[:,None],axis=0)
    base = float(min(risks[:2]))
    assert fit.fun <= min(risks)+1e-8
    records.append({"study_id":study_id,"patients":rows.strain_id.nunique(),"quantile":tau,
       "best_of_two_reference_losses":base,"best_observed_convex_loss":fit.fun,
       "maximum_observed_gain_over_best_reference":1-fit.fun/base,
       **{f"weight_{learner}":v for learner,v in zip(experts,fit.x[:k])}})
output = Path("results/scratch/clinical_followup_20261106")
output.mkdir(parents=True,exist_ok=True)
table = pd.DataFrame(records)
table.to_csv(output/"observed_convex_bound.csv",index=False)
print(table[["study_id","maximum_observed_gain_over_best_reference"]].to_string(index=False))
