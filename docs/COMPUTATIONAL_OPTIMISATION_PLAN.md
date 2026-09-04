# Computational optimisation plan (registry B=200)

Status: design analysis (2026-09-01), before any further server spend.
Goal: cut the final-registry wall time from ~10-12 days to ~3-4 days without
touching any frozen statistical object (estimator, Dantzig program, mu rule,
CV structure, grids, seeds, targets, thresholds).

## 1. Measured cost anatomy (per task, n=200, p=500; P05 instrumentation)

| component | share | notes |
|---|---|---|
| Dantzig mu-CV: 4 coords x 8 mu x 4 folds = 128 dense LPs (CLARABEL IPM) | ~75-80% | 2.9-4.3 s/solve; O(p^3) per solve, 128 factorisations per task |
| first stage (score-L1 + refit + reprofile) | ~15-20% | already FISTA-optimised; cluster-parallelisable |
| one-step + sandwich + bookkeeping | <5% | |

At p=1000/2000 the Dantzig term grows O(p^3) (measured pilot rho: P07 n=800
p=500 ~30-40 min/task; the p=1000/2000 cells are the registry's heavy tail).

## 2. Prioritised levers (implementation layer only)

### Lever 1: persistent OSQP model with warm-start chaining — **REJECTED (2026-09)**
- Original idea: each fold's CV solves share the constraint matrix A and
  objective P across all mu; only the box bounds change.
- Realisation: osqp's C-level API (osqpSetup/osqpUpdate/osqpWarmStart,
  exposed inside the installed osqp 1.0.0 namespace) keeps ONE KKT
  factorisation per fold; every subsequent solve costs O((n+m)^2) with warm
  starts.
- Adjudication: decision-INCOMPATIBLE at every parameterisation tried for
  p >= 200. On the frozen smallest-feasible-mu rule, CLARABEL returns
  mu = 0.0251 (residual 0.0251189, 1.27 s single solve) while OSQP
  (eps=1e-8, polishing, 100k iterations, 14.8 s) still reports infeasible and
  grid selection falls to mu = 0.0501. This is intrinsic (first-order ADMM
  cannot drive the tight L1 box constraints to the required residual
  tolerance), not a tuning artefact. Full evidence, parameter sweep, and
  testing blind spots are in `COMPUTATIONAL_OPTIMISATION_LEVER1_ADJUDICATION.md`.
- The claimed 0.8 s / 20-50x figures were not reproducible and are withdrawn.
  The OSQP path remains in the tree (default-off; env `JPDMEQR_DANTZIG_SOLVER`
  defaults to `cvxr`) only as instrumented evidence.

### Direction 1: parallelise the mu-CV CLARABEL grid — **IMPLEMENTED (2026-09)**
- The dominant term is 128 dense LPs per task (4 coords x 8 mu x 4 folds);
  every (fold, mu) solve is an independent pure function of its inputs.
  Solving the grid concurrently and aggregating in the original fold order is
  bit-identical to the serial reference (same mean/sd inputs, same argmin-tie
  rule, same infeasible-fold drop rule).
- Selection: option/env `jpDMEQR.cv_cores` / `JPDMEQR_CV_CORES`; default 1 =
  frozen serial path (never entered when <= 1). Linux: fork-based mclapply
  inside the task (keep `jobs` at 1 and use shards for task-level parallelism -
  no nested fork, per the PAPP-LD memory lesson). Windows: deterministic PSOCK
  with dynamic load balancing (`parLapplyLB`; CLARABEL wall time varies
  2.9-7.4 s per cell at p=500, so static round-robin would strand each fold's
  whole mu chain on one worker).
- Measured locally (8-core Windows desktop, p=500, 4 workers): serial 136.8 s
  vs parallel 49.7 s per coordinate-CV -> 2.75x; selected mu, min defect, and
  every candidate row bit-identical. p=120: 1.93x.
- Correctness gate: `tests/testthat/test-cv-parallel-parity.R` asserts
  bit-identity between serial and cv_cores=2/4 across coordinates.
- Registry effect at fixed core count: wall time of the dominant CV term
  divides by `cv_cores` (e.g. jobs=12 x cv_cores=4 on 48 logical cores) with
  no change to any statistical object. Independent further gains remain open:
  Lever 3 (server OpenBLAS 2-4x) and Lever 2 (1.2-1.5x).

### Lever 2: share the first-stage fit across methods within a replication
- PROFILE-DQR / POP-H / TRUE-SUPPORT re-run the identical calibrated first
  stage (L1 at h_est + post-refit at h_inf) three times per replication in
  current adapters. A shared-fit context (same fit object, deterministic)
  removes the duplication; only the precision/variance step differs.
- Expected effect: 1.2-1.5x on multi-method cells.
- Pending: dispatcher wiring (implementation only; results bit-identical).

### Lever 3: threaded BLAS on the server (free 2-4x on all dense ops)
- The Ubuntu reference BLAS is single-threaded; installing/libreplacing with
  OpenBLAS (apt libopenblas0 + update-alternatives) speeds every dense matmul
  (first stage Hessian assembly, CLARABEL dense ops).
- Expected effect: 2-4x across the whole pipeline; zero statistical change.
- Note: BLAS threading interacts with the process-count scheduling; tune the
  thread count (OPENBLAS_NUM_THREADS) to ~2-4 with 16-32 tasks in flight.

### Lever 4 (parked, documented for the paper): parametric path / ADMM
- O(p^2)-per-path homotopy for the Dantzig row; earlier PDHG experiment hit a
  degenerate saddle for loose mu (fixed-point analysis completed; epigraph
  form escape requires care). Parked - not needed to reach the 3-4 day target;
  a valid future-work paragraph for the manuscript's computational remarks.

## 3. Frozen invariants (must not change)

- Dantzig program min ||omega||_1 s.t. ||H omega - e_k||_inf <= mu;
- mu grid values and the smallest-feasible rule; the 2/4-fold CV structure and
  its loss; all seeds; all registry rows/config/checkpoint versions;
- first-stage formula, bandwidth h_est/h_inf, refit set rule (round 5);
- score/Hessian/sandwich formulas; gate thresholds and V1 semantics.

## 4. Costs

- Lever 1+2 implementation + parity + geometry re-validation: ~0.5-1 day
  (local; parallel with asset fetch if desired).
- One gate-revalidation cycle on the new commit (~1 h, geometry the long pole).
- No effect on the frozen registry contents (same B=200 rows).

## 5. Expected schedule after approval

1. implement levers 1+2 (+ OpenBLAS on the server at launch);
2. parity + strict geometry on the new commit;
3. brief server power-on: fetch assets AND install OpenBLAS;
4. freeze preflight (b200 config) -> launch registry at ~64 cores:
   **~3-4 days wall** (vs 10-12 without the plan);
5. application layer + manuscript backfill in parallel.

Sign-off requested from the study lead; theory visibility for the record
(solver substitution and intra-replication fit sharing are implementation
layer; no statistical object changes; the amendment file notes this plan).