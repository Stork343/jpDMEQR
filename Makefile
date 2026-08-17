.PHONY: install registry registry-all geometry geometry-strict pilot targets-dev targets-final micro-preflight freeze main summarise download-data build-data audit-data fit-data summarise-data manuscript

install:
	Rscript scripts/00_install_dependencies.R

registry:
	Rscript scripts/simulation/00_validate_registry.R --config=config/simulation_main.csv --type=main

registry-all:
	Rscript scripts/simulation/00_validate_registry.R --config=config/simulation_main.csv --type=main
	Rscript scripts/simulation/00_validate_registry.R --config=config/simulation_pilot.csv --type=pilot
	Rscript scripts/simulation/00_validate_registry.R --config=config/simulation_preflight.csv --type=preflight

geometry:
	Rscript scripts/simulation/00_validate_profile_geometry.R --strict=false

geometry-strict:
	Rscript scripts/simulation/00_validate_profile_geometry.R --strict=true --seed-count=20

pilot:
	Rscript scripts/simulation/01_run_pilot.R --config=config/simulation_pilot.csv --jobs=1

targets-dev:
	Rscript scripts/simulation/04_build_profile_targets.R --config=config/simulation_main.csv --development=true --n-population=10000 --repeats=2 --direction-repeats=2

targets-final:
	Rscript scripts/simulation/04_build_profile_targets.R --config=config/simulation_main.csv --n-population=100000 --repeats=4 --direction-repeats=4

micro-preflight:
	Rscript scripts/simulation/04_build_profile_targets.R --config=config/simulation_preflight.csv --development=true --n-population=5000 --repeats=2 --direction-repeats=2
	Rscript scripts/simulation/02_run_main.R --config=config/simulation_preflight.csv --development=true --allow-missing-benchmarks=true --jobs=1 --run-id=preflight_$$(date +%Y%m%dT%H%M%S)

freeze:
	Rscript scripts/simulation/05_preflight_freeze.R --config=config/simulation_main.csv

main:
	Rscript scripts/simulation/02_run_main.R --config=config/simulation_main.csv --jobs=1

summarise:
	Rscript scripts/simulation/03_summarise.R

download-data:
	Rscript scripts/application/00_download_GSE65391.R --processed=true --raw=false

build-data:
	Rscript scripts/application/01_build_GSE65391.R

audit-data:
	Rscript scripts/application/02_audit_GSE65391.R

fit-data:
	Rscript scripts/application/03_fit_GSE65391.R

summarise-data:
	Rscript scripts/application/04_summarise_GSE65391.R

manuscript:
	cd manuscript && bash compile.sh
