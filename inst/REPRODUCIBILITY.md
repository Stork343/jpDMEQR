# Reproducibility Checklist

This package ships with fixed workflows for simulation and real-data analyses used in the manuscript.

## Environment

- R version: 4.5.x or newer
- Required packages: `MASS`, `quantreg`
- Optional benchmark/real-data packages: `BGLR`, `lqmm`, `glmmLasso`, `GEOquery`, `Biobase`, `CVXR`

## Deterministic setup

- Global seed used in manuscript scripts: `set.seed(2026)`
- Simulation replication count for main tables: `B = 100`
- Six error scenarios: `S1`--`S6`

## Scripts

- Simulation: `inst/scripts/run_simulation_study.R`
- Real data: `inst/scripts/run_real_data_analysis.R`

## Core outputs

- Simulation raw and summary CSV files in `simulation_outputs/`
- Real-data combined table and diagnostics CSV files

## Verification items

1. `R CMD check` returns status `OK`.
2. Main simulation script finishes without runtime errors.
3. Output tables include estimation, selection, coverage, and runtime diagnostics.
4. Real-data script reports sample-size metadata and quantile-specific effect tables.
