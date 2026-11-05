import argparse
import itertools
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import xgboost as xgb
import yaml
from statsmodels.stats.weightstats import DescrStatsW

parser = argparse.ArgumentParser()
parser.add_argument("--config", default="config/application_charite_hb_20261105.yml")
parser.add_argument("--mode", choices=["raw", "offset"], required=True)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()
config = yaml.safe_load(Path(args.config).read_text())
root = Path(config["output_root"])
quantiles = np.asarray(config["quantiles"], dtype=float)


def read_export(key, part):
    directory = root / "exports_full_precision" / args.mode / key
    X = pd.read_csv(directory / f"{part}_X.csv")
    metadata = pd.read_csv(directory / f"{part}_meta.csv")
    scaling = json.loads((directory / "scaling.json").read_text())
    assert list(X.columns) == scaling["columns"]
    assert len(X) == len(metadata) and not metadata.sample_id.duplicated().any()
    assert np.isfinite(X.to_numpy()).all() and np.isfinite(metadata.y).all()
    assert np.isclose(metadata.weight.sum(), metadata.strain_id.nunique())
    return X, metadata, scaling


def matrix(X, metadata, reference=None):
    weights = len(metadata) * metadata.weight.to_numpy() / metadata.strain_id.nunique()
    return xgb.QuantileDMatrix(X, metadata.y, weight=weights, ref=reference,
                              max_bin=config["xgboost"]["max_bin"], nthread=1)


def parameters(setting, seed):
    return {**setting, "objective": "reg:quantileerror", "quantile_alpha": quantiles.tolist(),
            "tree_method": "hist", "device": "cpu", "multi_strategy": "one_output_per_tree",
            "subsample": config["xgboost"]["subsample"], "max_bin": config["xgboost"]["max_bin"],
            "nthread": 1, "seed": seed, "verbosity": 0}


def response_prediction(model, X, metadata, scaling, rounds):
    prediction = model.inplace_predict(X, iteration_range=(0, int(rounds)), strict_shape=True)
    assert prediction.shape == (len(metadata), len(quantiles)) and np.isfinite(prediction).all()
    return prediction * scaling["y_scale"] + scaling["y_center"] + metadata.offset.to_numpy()[:, None]


def calibrate(metadata, prediction):
    result = []
    for k, value in enumerate(quantiles):
        residual = metadata.raw_y.to_numpy() - prediction[:, k]
        calibration = float(DescrStatsW(residual, weights=metadata.weight).quantile([value], return_pandas=False)[0])
        below = np.average(residual < calibration - 1e-10, weights=metadata.weight)
        at_or_below = np.average(residual <= calibration + 1e-10, weights=metadata.weight)
        assert below <= value + 1e-8 and at_or_below >= value - 1e-8
        result.append(calibration)
    return np.asarray(result)


grid = [dict(zip(["max_depth", "learning_rate", "reg_lambda"], values)) for values in itertools.product(
    config["xgboost"]["max_depth"], config["xgboost"]["learning_rate"], config["xgboost"]["reg_lambda"])]

if args.check:
    X, metadata, scaling = read_export("inner_1", "fitting")
    V, held, held_scaling = read_export("inner_1", "validation")
    assert not set(metadata.strain_id) & set(held.strain_id)
    dtrain = matrix(X, metadata)
    dtest = matrix(V, held, reference=dtrain)
    metrics = {}
    model = xgb.train(parameters(grid[0], 20261105), dtrain, num_boost_round=20,
                      evals=[(dtest, "held")], evals_result=metrics, verbose_eval=False)
    native = model.predict(dtest, strict_shape=True)
    direct = model.inplace_predict(V, strict_shape=True)
    prediction_error = float(np.max(np.abs(native - direct)))
    assert prediction_error < 1e-6
    residual = held.y.to_numpy()[:, None] - native
    loss = np.maximum(quantiles * residual, (quantiles - 1) * residual)
    metric_error = abs(float(np.average(loss.mean(axis=1), weights=held.weight)) - metrics["held"]["quantile"][-1])
    assert metric_error < 1e-6
    path = root / "checks" / f"boost_{args.mode}.json"
    model.save_model(path)
    restored = xgb.Booster(model_file=path)
    reload_error = float(np.max(np.abs(restored.inplace_predict(V, strict_shape=True) - direct)))
    assert reload_error == 0
    calibrate(metadata, response_prediction(model, X, metadata, scaling, 20))
    output = Path(config["report_directory"]) / f"boost_{args.mode}_checks.json"
    output.write_text(json.dumps({"status": "passed", "prediction_error": prediction_error,
        "metric_error": metric_error, "reload_error": reload_error,
        "xgboost_version": xgb.__version__, "formal_validation_used": False}, indent=2) + "\n")
    print(output.read_text())
    raise SystemExit(0)

method = "XGB_RAW" if args.mode == "raw" else "XGB_OFFSET"
directory = root / "models" / f"{method}_all"
assert not directory.exists()
directory.mkdir(parents=True)
checks = json.loads((Path(config["report_directory"]) / f"boost_{args.mode}_checks.json").read_text())
assert checks["status"] == "passed"
settings = [dict(grid_index=g, rounds=r, calibration=mode) for g in range(len(grid))
            for r in config["xgboost"]["rounds"] for mode in config["xgboost"]["calibration"]]
losses = np.zeros((len(settings), len(quantiles)))
started = time.perf_counter()
(directory / "request.json").write_text(json.dumps({"method": method, "stage": "training", "config": config,
    "implementation_directory": str(Path(__file__).resolve().parent), "formal_validation_used": False}, indent=2) + "\n")
for fold in range(1, config["inner_folds"] + 1):
    X, metadata, scaling = read_export(f"inner_{fold}", "fitting")
    V, held, held_scaling = read_export(f"inner_{fold}", "validation")
    assert not set(metadata.strain_id) & set(held.strain_id)
    dtrain = matrix(X, metadata)
    for g, setting in enumerate(grid):
        print(f"{method} fold {fold} grid {g + 1}", flush=True)
        model = xgb.train(parameters(setting, 20261105 + fold), dtrain,
                          num_boost_round=max(config["xgboost"]["rounds"]), verbose_eval=False)
        model.save_model(directory / f"inner_{fold}_grid_{g}.json")
        predictions = []
        calibrations = []
        for rounds in config["xgboost"]["rounds"]:
            predicted = response_prediction(model, V, held, held_scaling, rounds)
            calibration = calibrate(metadata, response_prediction(model, X, metadata, scaling, rounds))
            predictions.append(predicted)
            calibrations.append(calibration)
            for s, item in enumerate(settings):
                if item["grid_index"] != g or item["rounds"] != rounds:
                    continue
                adjusted = predicted + (calibration if item["calibration"] == "training" else 0)
                residual = held.raw_y.to_numpy()[:, None] - adjusted
                values = np.maximum(quantiles * residual, (quantiles - 1) * residual)
                losses[s] += (values * held.weight.to_numpy()[:, None]).sum(axis=0)
        np.savez_compressed(directory / f"inner_{fold}_grid_{g}_predictions.npz",
            prediction=np.stack(predictions), calibration=np.stack(calibrations), sample_id=held.sample_id.to_numpy(dtype=str))

selected = losses.argmin(axis=0)
X, metadata, scaling = read_export("training", "fitting")
dtrain = matrix(X, metadata)
models = {}
calibrations = np.zeros(len(quantiles))
for g in sorted({settings[s]["grid_index"] for s in selected}):
    indices = [k for k, s in enumerate(selected) if settings[s]["grid_index"] == g]
    rounds = max(settings[selected[k]]["rounds"] for k in indices)
    model = xgb.train(parameters(grid[g], 20261115), dtrain, num_boost_round=rounds, verbose_eval=False)
    filename = f"final_grid_{g}.json"
    model.save_model(directory / filename)
    models[str(g)] = filename
    for k in indices:
        item = settings[selected[k]]
        if item["calibration"] == "training":
            prediction = response_prediction(model, X, metadata, scaling, item["rounds"])
            calibrations[k] = calibrate(metadata, prediction)[k]
fit = {"method": method, "model_type": "boost", "offset_model": args.mode == "offset",
       "quantiles": quantiles.tolist(), "models": models, "settings": settings, "grid": grid,
       "selected": selected.tolist(), "calibration": calibrations.tolist(), "loss": losses.tolist(),
       "columns": list(X.columns), "scaling": scaling, "xgboost_version": xgb.__version__}
(directory / "fit.json").write_text(json.dumps(fit, indent=2) + "\n")
(directory / "status.json").write_text(json.dumps({"status": "completed", "stage": "training",
    "runtime_sec": time.perf_counter() - started, "formal_validation_used": False, "manuscript_adoption": False}, indent=2) + "\n")
print(f"{method} training completed", flush=True)
