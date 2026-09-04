# jpDMEQR simulation-freeze ordinary-file patch

Target repository: `Stork343/jpDMEQR`

Target branch: `gpt/simulation-freeze-resolution`

This bundle contains ordinary UTF-8 source/configuration files. It is not a segmented archive and does not depend on the deleted `tools/simulation-freeze-patch/part*` files.

## Registry sizes

- `config/simulation_main.csv`: 79 physical lines (1 header + 78 registered scenarios)
- `config/simulation_pilot.csv`: 5 physical lines (1 header + P01--P04)
- `config/simulation_preflight.csv`: 11 physical lines (1 header + PF01--PF10)

## Files to replace or add

- `config/simulation_main.csv`
- `config/simulation_pilot.csv`
- `config/simulation_preflight.csv`
- `R/v2_simulation_v2_part01.R`
- `R/v2_simulation_v2_part02.R`
- `R/v2_simulation_v2_part03.R`
- `R/v2_benchmark_adapters_v2_part03.R`
- `scripts/00_install_dependencies.R`
- `scripts/simulation/00_validate_profile_geometry.R`
- `scripts/simulation/00_validate_registry.R`
- `scripts/simulation/01_run_pilot.R`
- `scripts/simulation/02_run_main.R`
- `scripts/simulation/03_summarise.R`
- `scripts/simulation/04_build_profile_targets.R`
- `scripts/simulation/05_preflight_freeze.R`
- `scripts/simulation/_run_registry_helpers_part01.R`
- `scripts/simulation/_run_registry_helpers_part02.R`
- `scripts/simulation/_run_registry_helpers_part03.R`
- `scripts/simulation/_run_registry_helpers_part04.R`
- `Makefile`

`SHA256SUMS` records content hashes for all patch files.

## Validation completed before packaging

- CSV dimensions, unique experiment IDs, and frozen scenario coverage checks passed.
- `simulation_main.csv` contains the complete Module A rate grid, error/quantile additions, signal-strength cells, copula and cluster-size imbalance cells, profile-target misspecification cells, all screening routes, the full nuisance-penalty grid, and computation-scaling cells.
- Static lexical checks found balanced R delimiters and strings in every supplied R file.

## Validation still required after repository application

Run the repository's R parse/unit tests and GitHub Actions. The construction environment used to package these files did not contain R, so this bundle does not claim a completed R-runtime validation.
