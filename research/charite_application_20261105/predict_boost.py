import json
from pathlib import Path

import numpy as np
import pandas as pd
import xgboost as xgb
import yaml

config = yaml.safe_load(Path("config/application_charite_hb_20261105.yml").read_text())
root = Path(config["output_root"])
release = json.loads((root / "validation_release.json").read_text())
assert release["status"] == "frozen"
records = []
checks = []
for mode in ["raw", "offset"]:
    method = "XGB_RAW" if mode == "raw" else "XGB_OFFSET"
    directory = root / "models" / f"{method}_all"
    fit = json.loads((directory / "fit.json").read_text())
    source = root / "validation" / "boost_features" / mode
    X = pd.read_csv(source / "X.csv")
    metadata = pd.read_csv(source / "metadata.csv")
    scaling = json.loads((source / "scaling.json").read_text())
    assert "y" not in metadata.columns and "raw_y" not in metadata.columns
    assert list(X.columns) == fit["columns"] == scaling["columns"]
    assert np.isclose(scaling["y_scale"], fit["scaling"]["y_scale"], rtol=0, atol=1e-12)
    assert np.isclose(scaling["y_center"], fit["scaling"]["y_center"], rtol=0, atol=1e-12)
    assert fit["xgboost_version"] == xgb.__version__
    matrix = xgb.DMatrix(X, nthread=1)
    models = {g: xgb.Booster(model_file=directory / filename) for g, filename in fit["models"].items()}
    for k, tau in enumerate(fit["quantiles"]):
        setting = fit["settings"][fit["selected"][k]]
        model = models[str(setting["grid_index"])]
        raw = model.inplace_predict(X, iteration_range=(0, setting["rounds"]), strict_shape=True)[:, k]
        reference = model.predict(matrix, iteration_range=(0, setting["rounds"]), strict_shape=True)[:, k]
        error = float(np.max(np.abs(raw - reference)))
        assert error < 1e-6
        prediction = raw * scaling["y_scale"] + scaling["y_center"] + metadata.offset + fit["calibration"][k]
        assert np.isfinite(prediction).all()
        records.append(pd.DataFrame({"method": method, "sample_id": metadata.sample_id,
            "strain_id": metadata.strain_id, "tau": tau, "prediction": prediction}))
        checks.append({"method": method, "tau": tau, "prediction_interface_error": error})
pd.concat(records, ignore_index=True).to_csv(root / "validation" / "predictions_boost.csv", index=False)
pd.DataFrame(checks).to_csv(Path(config["report_directory"]) / "boost_prediction_checks.csv", index=False)
