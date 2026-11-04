import json
import subprocess
from pathlib import Path

import yaml

root = Path("data/derived/yeast_1011_20261104")
logs = Path("results/logs/yeast_1011_20261104")
logs.mkdir(parents=True, exist_ok=True)
config = yaml.safe_load(Path("config/discovery_yeast_1011_20261104.yml").read_text())
groups = json.loads((root / "training_sets.json").read_text())
plink = "results/scratch/dependencies/plink2/plink2"
source = "data/raw/yeast_1011_20261104/gwas/1011GWAS_matrix"
for name in groups:
    prefix = root / f"{name}_pruning"
    panel = root / f"{name}_panel"
    assert not panel.with_suffix(".bed").exists()
    with (logs / f"{name}_pruning.log").open("w") as log:
        subprocess.run([
            plink, "--bfile", source, "--chr-set", "17", "no-x", "no-y", "no-xy", "no-mt", "--chr", "1-16",
            "--keep", str(root / f"{name}.keep"), "--maf", str(config["marker_filter"]["training_maf_minimum"]),
            "--geno", str(config["marker_filter"]["training_missing_fraction_maximum"]),
            "--indep-pairwise", f'{config["marker_filter"]["ld_window_kb"]}kb', "1", str(config["marker_filter"]["ld_r_squared"]),
            "--threads", "2", "--memory", "4096", "--out", str(prefix),
        ], stdout=log, stderr=subprocess.STDOUT, check=True)
    with (logs / f"{name}_panel.log").open("w") as log:
        subprocess.run([
            plink, "--bfile", source, "--chr-set", "17", "no-x", "no-y", "no-xy", "no-mt", "--chr", "1-16",
            "--keep", str(root / "eligible.keep"), "--extract", str(prefix) + ".prune.in", "--make-bed",
            "--threads", "2", "--memory", "4096", "--out", str(panel),
        ], stdout=log, stderr=subprocess.STDOUT, check=True)
    print(name, "completed", flush=True)
