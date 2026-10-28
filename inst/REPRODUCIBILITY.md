# Reproducibility protocol for the SJS reconstruction

The active protocol is defined by `docs/METHOD_SPECIFICATION.md`,
`docs/SIMULATION_PROTOCOL.md`, `docs/SIMULATION_FREEZE_DECISIONS.md`,
`docs/RESULTS_CONTRACT.md`, `docs/REGISTRY_AMENDMENT_B200.md`, the registries
in `config/`, and, for the real-data layer,
`docs/EMPIRICAL_PROTOCOL_GSE65391.md`. Manuscript numbers come only from the
frozen B=200 registry `config/simulation_main_b200.csv` executed under a
matching freeze manifest. The earlier B=100 six-scenario benchmark is
historical.

The frozen pipeline (registries, runners, governance files) lives in the
repository checkout; the installed binary package ships interfaces only, so
reproduction always runs from a checkout. Every tier below runs on an ordinary
multi-core desktop with 32 GB of memory. No server is required for any tier
except the full manuscript-scale registry, whose schedule is documented for
both a server and a sharded desktop.

Prerequisites: R 4.6.1 and the dependencies installed by
`scripts/00_install_dependencies.R`. The `dual_warm` Dantzig path additionally
needs the optional `highs` package from the same script.

## Solver configuration switches

Six environment variables control the implementation layer. Every default
reproduces the frozen reference path exactly; none of them touches a frozen
statistical object (estimator, Dantzig program, mu rule, CV structure, grids,
seeds, targets, thresholds).

| switch | default | effect |
|---|---|---|
| `JPDMEQR_DANTZIG_SOLVER` | unset (frozen CVXR path) | `dual_warm` solves the dual program on a persistent HiGHS simplex model with the basis carried along the mu grid: 4.05-6.00x on the dominant mu-CV term at p=500, 17.07x serial and 30.84x over four folds at p=1000, above 36x at p=2000. Same selected mu, same fold counts, inference columns agree below 1.1e-08 relative |
| `JPDMEQR_CV_CORES` | 1 | fold-parallel mu-CV; aggregates in fold order, bit-identical to serial (`test-cv-parallel-parity.R`) |
| `JPDMEQR_CLUSTER_SOLVER` | bfgs | `newton` uses the analytic damped Newton nuisance solve: 1.47x on one stage-1 fit at A02, 1.30x at A08 |
| `JPDMEQR_ASSEMBLY` | sequential | `batch` assembles Hessians with whole-matrix products: 1.19x at A02 up to 2.21x at p=1000 |
| `JPDMEQR_HESSIAN_REFRESH` | 5 | first-stage curvature refresh interval; 20 gives 1.03-1.28x with fits inside the frozen 1e-3 acceptance. The round-5 refit, the jackknife, the population target assets and the geometry audit keep the frozen 5 |
| `JPDMEQR_SHARED_FIRST_STAGE` | off | PROFILE-DQR, PROFILE-DQR-POP-H and POSTREFIT-EXACT-H reuse one calibrated first-stage fit per replication; bit-identical to independent fits (`test-shared-first-stage.R`). Saves one stage-1 fit per additional sharing method per replication |

Measured whole-task effects of the stack: one PROFILE-DQR replication falls
from 569.6-651.6 s to 98.0-100.0 s at A02 (p=500), from 3825.0 s to 220.8 s at
A05 (p=1000) and from more than 7.7 h to 729.2 s at A08 (p=2000) on the
development server; details and adjudications are in
`docs/SERVER_RUNBOOK.md` and `docs/COMPUTATIONAL_OPTIMISATION_PLAN.md`.

Recommended local stack. Set these variables once per shell and use the same
set for the geometry gate and for every registry run:

PowerShell:

```powershell
$env:JPDMEQR_DANTZIG_SOLVER = 'dual_warm'
$env:JPDMEQR_CLUSTER_SOLVER = 'newton'
$env:JPDMEQR_ASSEMBLY = 'batch'
$env:JPDMEQR_HESSIAN_REFRESH = '20'
$env:JPDMEQR_SHARED_FIRST_STAGE = '1'
$env:JPDMEQR_CV_CORES = '4'
```

bash:

```bash
export JPDMEQR_DANTZIG_SOLVER=dual_warm JPDMEQR_CLUSTER_SOLVER=newton
export JPDMEQR_ASSEMBLY=batch JPDMEQR_HESSIAN_REFRESH=20
export JPDMEQR_SHARED_FIRST_STAGE=1 JPDMEQR_CV_CORES=4
```

`solver_configuration_v2()` records all six values in the geometry manifest,
the freeze manifest and every run directory (`solver_configuration.json`). A
gate passed under one configuration does not authorise another: changing any
variable between the gate and a run requires re-running the geometry gate in
the new configuration, and the entry scripts stop with that instruction. At
p=2000, `JPDMEQR_CV_CORES=2` is the safe ceiling on a 32 GB machine
(`docs/SERVER_RUNBOOK.md`).

## Tier one unit tests and the strict geometry gate

From the repository root:

```
Rscript -e "pkgload::load_all('.'); testthat::test_dir('tests/testthat')"
Rscript scripts/simulation/00_validate_profile_geometry.R --strict=true
```

The suite covers the profile gradient, Hessian and Schur identities, the
Dantzig row program and the mu selection rule, solver parity (fold-parallel
against serial, dual_warm against the frozen path), the shared-first-stage
cache, every benchmark adapter and the configuration-manifest logic. It
measured 9.0 minutes on an 8-core Windows desktop under the recommended stack
(longer under the frozen defaults, where the CVXR-path parity tests dominate).

The strict geometry gate runs 760 finite-difference cases over 20 seeds across
all target quantiles, q=1/2, h multipliers, balanced and unbalanced clusters
and non-diagonal Lambda, plus step-size sensitivity, the mandatory Dantzig
validation and the assembly-mode identity comparison. It measured 3.6 minutes
under the recommended stack on the same desktop (the frozen-default design
estimate is about an hour). It writes `results/preflight/geometry/manifest.json`
bound to the current commit and the current configuration; every registry entry
point re-checks that manifest before running anything.

## Tier two micro preflight

```
Rscript scripts/simulation/02_run_main.R --config=config/simulation_preflight.csv --run-id=micro_preflight_local --development
```

PF01-PF10, B=2, 24-30 clusters, p=30-60: every DGP branch, every screening
route, the target loaders, adapter dispatch and the output schemas. It
measured 5.1 minutes on the same desktop under the recommended stack, writes
the micro-preflight manifest and stops on any runner failure or any method row
that is not ok.

## Tier three quick reproduction cells

```
Rscript scripts/simulation/09_quick_repro.R
```

`Rscript inst/scripts/run_simulation_study.R` forwards to the same entry from
any directory inside the checkout. Both accept `--run-id`, `--jobs`,
`--max-reps` and `--experiment-id`. On POSIX `--jobs` parallelises the
replications by fork; on Windows the runner keeps the replications serial
inside one process and task-level parallelism goes through sharded launches
(`docs/DESKTOP_RUNNER.md`), while the in-task fold parallelism from
`JPDMEQR_CV_CORES` is active on both platforms. The entry enforces the
governance files, the registry contract and the geometry gate, then runs
`config/simulation_quick_repro.csv`:

| cell | design | B | methods |
|---|---|---|---|
| QR01 | A02 design, n=200, p=500, s=5, tau=0.5, structural target | 5 | PROFILE-DQR, POSTREFIT-EXACT-H, POOLED-QR-LASSO, SQR-DEBIASED-IID, PROFILE-DQR-TRUE-SUPPORT |
| QR02 | A05 design, n=200, p=1000, structural target | 2 | PROFILE-DQR, POSTREFIT-EXACT-H |

This is a demonstration tier: it verifies that the frozen pipeline runs
end-to-end on local hardware and lets a reviewer measure the per-cell cost.
Its outputs never enter the manuscript.

Measured on an 8-core Windows desktop (R 4.6.1), jobs=1, cv_cores=4,
recommended stack, per replication: PROFILE-DQR 149 s at p=500 and 443 s at
p=1000; SQR-DEBIASED-IID 400 s; PROFILE-DQR-TRUE-SUPPORT 29 s; POOLED-QR-LASSO
2.5 s; POSTREFIT-EXACT-H below 1 s with the shared first stage on, because it
reuses the PROFILE-DQR fit object. The complete quick reproduction takes about
50 minutes for QR01 and about 15 minutes for QR02 on this machine. Under the
frozen defaults a PROFILE-DQR replication costs 5.8-6.5x more at p=500 and
17.3x more at p=1000, so the same tier extends to several hours.

Summaries: `Rscript scripts/simulation/03_summarise.R --run-id=<run_id>`.

## Full manuscript scale reproduction

The frozen registry `config/simulation_main_b200.csv` (78 rows, B=200) runs
only through the complete validation sequence in `AGENTS.md` section 10:
benchmark acceptance manifests on the current commit, final population target
assets (`scripts/simulation/04_build_profile_targets.R`, at least 100,000
clusters and four independent repeats), the pilot gate decision, and the freeze
preflight (`scripts/simulation/05_preflight_freeze.R`) under a manifest that
matches the current commit and registry checksum. The final execution is
`scripts/simulation/02_run_main.R` without any development flag. Sharded launch
pattern: `scripts/simulation/launch_sharded_pilot.R`; desktop porting
instructions: `docs/DESKTOP_RUNNER.md`. Server schedule with the switch stack:
about 3-4 days on 48-64 cores. A sharded desktop can reproduce individual
cells: one PROFILE-DQR replication of a p=500 cell costs 98-100 s under the
stack (about 5.5 core-hours for B=200), and the comparators in a cell add
their own per-replication cost on top (measured in `docs/SERVER_RUNBOOK.md`
and in the tier-three table above).

## Determinism

The base seed is 20260817. Each replication receives a deterministic integer
seed derived from the base seed, experiment identifier, replicate and role.
Parallel scheduling never changes generated data or results: fold-parallel
mu-CV aggregates in fold order and is bit-identical to serial, and the
shared-first-stage cache compares full request lists with `identical()` and
hands out the same deterministic fit object. Outputs are sorted deterministically
by experiment, replicate, method and coordinate.

## No calibration rule

Tuning may depend on training data, theory-motivated grids, cluster-level
cross-validation and numerical feasibility. It may not depend on the true
coefficient, empirical coverage, empirical bias or held-out test outcomes.
Failed fits remain in the denominator and failure rates are reported.

## Legacy interfaces

The pre-reconstruction exports (`Simulation_Study`,
`export_simulation_outputs`, `run_real_data_analysis`) remain installed for
continuity; the authoritative simulation pipeline is the v2 module set
(`R/v2_*.R`) with `scripts/simulation/`. `inst/scripts/run_simulation_study.R`
forwards to the v2 quick reproduction entry. Empirical work follows
`docs/EMPIRICAL_PROTOCOL_GSE65391.md` and `scripts/application/`: download,
build and audit are authorised, while confirmatory fitting waits for the
empirical gates in `docs/EMPIRICAL_FREEZE_DECISIONS.md`.
