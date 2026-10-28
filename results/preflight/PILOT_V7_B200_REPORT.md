# Pilot v7 gate decision

## Run identification

| Item | Value |
| --- | --- |
| run_id | pilot_v7_B200 (shards pilot_v7_B200_s1 to s4) |
| Implementation commit | 0b4023b1f3f472709283d4aa563f88847b9de82b |
| Branch | gpt/simulation-freeze-resolution |
| Registry checksum (config_sha256) | 915760ced523b24833888b406c5d8cb9cb169ed8fdf91214a58803f6ae7390e9 |
| Cells and replications | P01 to P07, B = 200 each, 1400 replications total |
| Execution window | started 2026-09-30 23:57:59 local; shards completed 2026-10-03 08:27 to 09:19 local |
| Execution style | 4 independent shards of 350 tasks, Windows in-process serial runner (jobs = 1), JPDMEQR_CV_CORES = 2 fold workers per shard |
| Hardware | HOUJIAN, Windows 10 x64, 8 logical cores, R 4.6.1 ucrt |
| Solver stack | dantzig_solver dual_warm, cluster_solver newton, assembly batch, cv_workers 2, hessian_refresh 5, shared_first_stage FALSE |

The current strict geometry manifest is bound to the same commit and the same solver stack, so the pilot ran on the geometry-certified configuration.

## Decision

**PASS** under the V1 mechanism criteria of the round-6 amendment, where the SE and coverage bands are recorded as visible diagnostics and carry no binary veto. The machine-readable manifest is `results/preflight/pilot_gate.json` (evaluated 2026-10-03 11:13:01 UTC, pass true, commit 0b4023b).

No quantity was tuned to obtain this decision. Thresholds were not touched, and the completed pilot authorises no post hoc tuning of h, lambda_beta, Lambda, mu, screening options, standard-error multipliers, or acceptance bands.

## Record completeness

| Table | Rows |
| --- | --- |
| replication_metrics | 8200 (8195 ok, 4 warning, 1 failed) |
| coordinate_metrics | 26396 |
| theory_diagnostics | 8200 |
| screening_records | 1400 |
| runner_failures | 0 |
| not_implemented rows | 0 |

The single failed replication is P05 replicate 99 of PROFILE-DQR, failure stage precision_solver, runtime 7520 s. The row is retained in the raw tables with its failure stage; it is neither dropped nor zero-encoded. Merged row counts equal the sum over the four shards for every table.

## Mechanism criteria for PROFILE-DQR

| Criterion | Threshold | Observed | Verdict |
| --- | --- | --- | --- |
| Success fraction per cell P01 to P06 | >= 0.98 | 1.000 for P01 to P04, P06, P07; 0.995 for P05 | pass |
| Profile identity error | median < 1e-8, max < 1e-6 | median 5.0e-16, max 1.94e-15 (1399 finite rows) | pass |
| Dantzig feasibility | >= 0.95 of coordinate rows | 1.000000 over 5596 rows, zero infeasible rows | pass |
| Runner failures | 0 | 0 | pass |
| Missing methods | none | none | pass |

## Diagnostics recorded without calibration

Small-n band diagnostics. At P05 (n = 200) the PROFILE-DQR standard-error to empirical-SD ratios are 0.27 to 0.59 with coverage 0.37 to 0.75; coordinate x00002 (active) shows bias +0.10 against an empirical SD of 0.13. At P06 (n = 400) the ratios improve to 0.66 to 0.78 with coverage 0.79 to 0.89. The monotone improvement with n is consistent with the round-6 variance-ladder findings that the plug-in sandwich variance understates sampling variability at intentionally small pilot n, which is the documented reason these bands are diagnostics.

Inverse-defect diagnostics. D_k median 5.80, q90 13.61. A_k median 0.47, q90 1.84. D_k is a Holder upper bound and may be loose; no hard threshold applies.

Numerical geometry diagnostics. Identity errors are at machine precision in every cell. Max nuisance gradient over all rows is 1.64e-8. Mean Hessian max error is nonzero at P01 (1.10e-1), P05 (7.93e-2), P06 (5.95e-2) with matching precision-row errors 4.80, 4.52, 3.64; it is exactly zero at P02 to P04 and P07.

Per-replication mean runtime in seconds for PROFILE-DQR: P01 25.6, P02 36.2, P03 34.3, P04 31.3, P05 357.6, P06 387.5, P07 884.5. These values feed the wall-time planning of the final run.

## Paired pinball differences against POOLED-QR-LASSO

| Method | Pairs | Mean difference | MCSE |
| --- | --- | --- | --- |
| PROFILE-DQR | 1199 | -0.174 | 0.0022 |
| SQR-DEBIASED-IID | 800 | -0.172 | 0.0022 |
| BIAS-ADJ-LQMM | 200 | -0.179 | 0.0035 |
| DOUBLE-PEN-QLMM | 400 | -0.076 | 0.0027 |
| LQMM | 400 | -0.242 | 0.0042 |
| QGEE-SCAD | 400 | -0.140 | 0.0049 |
| QIF-SEE | 200 | -0.063 | 0.0031 |
| POSTREFIT-EXACT-H | 600 | -0.178 | 0.0025 |
| PROFILE-DQR-POP-H | 600 | -0.178 | 0.0025 |
| PROFILE-DQR-TRUE-SUPPORT (oracle) | 1200 | -0.258 | 0.0026 |
| PROFILE-DQR-TRUE-NUISANCE (oracle) | 400 | -0.200 | 0.0047 |

Negative values mean lower pinball loss than the pooled baseline. Oracle-labelled rows are visible diagnostics and are excluded from any practical superiority claim.

## Stage status after this decision

Stage 6 (pilot) is complete: raw results, failure records, coordinate metrics, MCSE summaries, theory diagnostics, and this written gate decision exist.

Prerequisites that remain open before the Stage 7 freeze preflight:

1. Final-scale profile_mc target assets for C06, C07, C14, C15, C16, D01 to D06 are under construction (serial builder started 2026-10-03 09:25 local on the fixed builder script).
2. Population-direction assets for A01 to A12 at 100000 clusters and four repeats completed on 2026-10-02 and 2026-10-03 at commit 0b4023b.
3. PF04 to PF06 asset refresh at the current commit is queued for relaunch after the overnight process loss.

The pilot pass alone does not authorise the final run. Stage 8 starts only after a freeze manifest matches the current commit and registry checksum.
