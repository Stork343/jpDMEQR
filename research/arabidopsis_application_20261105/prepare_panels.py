import json
import subprocess
from pathlib import Path

import yaml

root = Path("data/derived/arabidopsis_20261105")
logs = Path("results/logs/arabidopsis_20261105")
config = yaml.safe_load(Path("config/discovery_arabidopsis_20261105.yml").read_text())
panels = root / config["genotype_panel_subdirectory"]
panels.mkdir(parents=True, exist_ok=True)
logs = logs / config["genotype_panel_subdirectory"]
logs.mkdir(parents=True, exist_ok=True)
groups = json.loads((root / "training_sets.json").read_text())
plink = "results/scratch/dependencies/plink2/plink2"
for name in groups:
    prefix = panels / f"{name}_pruning"
    panel = panels / f"{name}_panel"
    assert not panel.with_suffix(".bed").exists()
    with (logs / f"{name}_pruning.log").open("w") as log:
        subprocess.run([
            plink, "--bfile", str(root / "full_genotypes"), "--chr-set", "5", "no-x", "no-y", "no-xy", "no-mt",
            "--keep", str(root / f"{name}.keep"), "--maf", str(config["marker_filter"]["training_maf_minimum"]),
            "--indep-pairwise", f'{config["marker_filter"]["ld_window_kb"]}kb', "1", str(config["marker_filter"]["ld_r_squared"]),
            "--threads", "2", "--memory", "4096", "--out", str(prefix),
        ], stdout=log, stderr=subprocess.STDOUT, check=True)
    with (logs / f"{name}_panel.log").open("w") as log:
        subprocess.run([
            plink, "--bfile", str(root / "full_genotypes"), "--chr-set", "5", "no-x", "no-y", "no-xy", "no-mt",
            "--extract", str(prefix) + ".prune.in", "--make-bed", "--threads", "2", "--memory", "4096",
            "--out", str(panel),
        ], stdout=log, stderr=subprocess.STDOUT, check=True)
    print(name, "completed", flush=True)
