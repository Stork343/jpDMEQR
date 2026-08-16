# Reproducibility protocol: SJS reconstruction

The earlier `B=100` six-scenario benchmark is historical and is not evidence for the reconstructed estimator. The active protocol is defined by:

- `docs/METHOD_SPECIFICATION.md`
- `docs/SIMULATION_PROTOCOL.md`
- `config/simulation_pilot.csv`
- `config/simulation_main.csv`
- `docs/EMPIRICAL_PROTOCOL_GSE65391.md`
- `docs/RESULTS_CONTRACT.md`

## Required order

1. Install dependencies.
2. Pass all profile-geometry tests.
3. Run the pilot registry (`B=200`).
4. Review pilot acceptance diagnostics.
5. Freeze implementation and configurations.
6. Run core inference scenarios (`B=1000`) and secondary scenarios (`B=500`).
7. Download and audit GSE65391.
8. Freeze subject folds and preprocessing.
9. Run the empirical analysis.
10. Generate all manuscript tables and figures from raw outputs.

## Determinism

The base seed is `20260817`. Each replication receives a deterministic integer seed derived from the base seed, experiment identifier and replication number. Parallel scheduling must not change generated data or results.

## No-calibration rule

Tuning may depend on training data, theory-motivated grids, cluster-level cross-validation and numerical feasibility. It may not depend on the true coefficient, empirical coverage, empirical bias or held-out test outcomes. Failed fits remain in the denominator and failure rates are reported.
