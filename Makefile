.PHONY: install geometry pilot targets main summarise download-data build-data audit-data fit-data summarise-data manuscript

install:
	Rscript scripts/00_install_dependencies.R

geometry:
	Rscript scripts/simulation/00_validate_profile_geometry.R

pilot:
	Rscript scripts/simulation/01_run_pilot.R --config=config/simulation_pilot.csv --jobs=1

targets:
	Rscript scripts/simulation/04_build_profile_targets.R --config=config/simulation_main.csv --n-population=10000 --repeats=2

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
