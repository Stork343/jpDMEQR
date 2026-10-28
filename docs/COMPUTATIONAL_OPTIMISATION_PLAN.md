# Computational optimisation plan (registry B=200)

Status: all levers adjudicated. Direction 1 and Levers 2, 5, 7 and 8 are
implemented as opt-in switches (the defaults keep the frozen reference path);
Lever 1 is rejected, Lever 6 is adjudicated inadvisable, Lever 3 waits on the
server OpenBLAS install and Lever 4 stays parked. Design analysis 2026-09-01,
before any further server spend.
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
  no change to any statistical object. Further independent gains: Lever 3
  (server OpenBLAS 2-4x) remains open; Lever 2 is implemented below.

### Lever 2: share the first-stage fit across methods within a replication - IMPLEMENTED, opt-in
- Switch `JPDMEQR_SHARED_FIRST_STAGE` (option `jpDMEQR.shared_first_stage`),
  default off = frozen reference path. When on, the methods that request the
  identical calibrated first stage on the full design within one replication -
  PROFILE-DQR, PROFILE-DQR-POP-H and POSTREFIT-EXACT-H - reuse one fit object
  from a replication-scoped cache environment held by the registry runner.
  The calibrated first stage is deterministic (no RNG on its fit path), so the
  shared result is bit-identical to independent fits; only `runtime_sec` moves.
- Requests that differ by design never share: TRUE-SUPPORT refits on the
  oracle sub-design with `penalty_factor = 0`, TRUE-NUISANCE modifies the
  response and nuisance design, SQR-DEBIASED-IID skips calibration and SPLIT
  uses disjoint cluster subsets. The cache compares the full request list with
  `identical()`, so any difference computes a fresh fit.
- Saving: one calibrated stage-1 fit per additional sharing method per
  replication. Measured stage-1 costs under the switch stack: 40.31 s at A02
  and 417.97 s at A08 (`p=2000`); the frozen-default fits cost more
  (60.45 s and 544.03 s at the same cells).
- Governance: the switch is recorded by `solver_configuration_v2()` and bound
  into the geometry manifest, so adopting it requires a strict geometry gate
  run under the same configuration, like every other knob.
- Correctness gate: `tests/testthat/test-shared-first-stage.R` asserts
  bit-identical adapter outputs with and without sharing for all three
  methods, a single cache entry for identical requests, and dispatcher gating
  (the cache is injected only while the switch is on).

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

### Lever 5: batch assembly of the profile Hessian - IMPLEMENTED, opt-in
- `assembly = "batch"` (option `jpDMEQR.assembly`, environment variable
  `JPDMEQR_ASSEMBLY`) keeps w_i, X_i'W_iZ_i and A_i per cluster and assembles
  the full Hessian, the Schur correction and the four fold Hessians with one
  large `crossprod` per sum. The default `sequential` is unchanged and stays
  bit-identical.
- Measured effect: 2.21x at p=1000/n=100, 1.96x at p=2000/n=200, 1.72x at
  p=2000/n=400 on the component call; relative difference 1e-15; peak memory
  +6 percent at p=2000.
- Measured and rejected in the same examination: hoisting X_i'beta and
  memoising Z_i'gamma. At m_i=5 those products cost microseconds against
  hundreds of microseconds of per-cluster interpreter work, so the change is
  neutral in time and adds branches to every evaluation.
- See `COMPUTATIONAL_OPTIMISATION_LEVER5_REPORT.md`.

### Lever 6: solve the compact Dantzig LP with HiGHS - ADJUDICATED INADVISABLE
- `JPDMEQR_DANTZIG_SOLVER=highs` (option `jpDMEQR.dantzig_solver`) hands the
  identical compact form C of the frozen LP to HiGHS instead of CLARABEL. The
  A matrix, the objective, the row bounds, the mu grid and the smallest-feasible
  rule are unchanged; only the solver changes.
- On well-conditioned constructions of the same LP the swap is a clear win:
  2.37x at p=300 and 3.12x at p=600 on one row over the frozen grid, and
  2.02x/2.70x at the mu-CV level, with the same selected mu, the same fold
  counts and defect tables agreeing to 2.3e-08 and 4.3e-07.
- On the real profile Hessian the win reverses. At the registered cell A02
  (n=200, p=500) one Dantzig row costs the compact CLARABEL path 15.8 s and
  15.0 s, against 181.3 s and 181.6 s for HiGHS, which hit a 60 s per-solve
  limit on three of the four candidates. The real Hessian is numerically
  singular, smallest singular value 2.589e-19 against a largest of 4.153e+00,
  while the construction used for the timing comparison has condition 6.04.
  Both paths still select the same mu (1.079702e-01).
- The construction of the constraint matrix accounts for 0.6 to 1.2 percent of
  a solve; the remaining cost is solver iterations, and on this geometry the
  interior-point code absorbs the numerical trouble better. The `start` argument
  of `highs_solve`, which passes a primal point to a one-shot call, gave no gain
  on the well-conditioned cells.
- The path stays in the code behind its switch as instrumented evidence and is
  not used by any run. See `COMPUTATIONAL_OPTIMISATION_HIGHS_PATH_REPORT.md`.

### Lever 7: dual form on a persistent simplex model - IMPLEMENTED (2026-09)
- `JPDMEQR_DANTZIG_SOLVER=dual_warm` solves the dual of the frozen row program,
  `min e_k'w + mu||w||_1` subject to `||H'w||_inf <= 1`, on one persistent HiGHS
  simplex model per fold and coordinate. The constraint matrix
  `rbind(cbind(H',-H'), cbind(-H',H'))` and the right-hand side depend on `H`
  alone; only the cost vector moves with `mu` and with the coordinate.
- The basis is carried along the descending `mu` grid, so the eight candidates
  are reoptimisations of one model: 0, 1, 0 and 12 simplex iterations on the four
  large candidates of A02 fold 1, 1085 on the candidate that binds. A candidate
  that returns an unbounded certificate settles every smaller candidate without
  a solve, and those settlements are recorded in the attempts table.
- Measured at the registered cell A02 on the real fold Hessians: the mu-CV call
  for coordinate 1 falls from 115.55 s to 28.53 s and for coordinate 3 from
  113.49 s to 18.92 s, a factor of 4.05 and 6.00. Both select the same `mu`
  (`4.3188077227e-01`), the same five eligible candidates and the same fold
  counts; the defect tables agree to 7.7e-09 and 6.2e-07.
- The four folds share no state, so they are the parallel unit: with
  `JPDMEQR_CV_CORES=4` the same two calls take 8.81 s and 5.75 s, against 30.15 s
  and 23.13 s for their own serial runs, a factor of 3.42 and 4.02. Each fold is
  a pure function of its training matrix, so the parallel table equals the serial
  table bit for bit. Together with the basis chain the call falls from 115.55 s to
  8.81 s and from 113.49 s to 5.75 s.
- Sharing one model across the four coordinates of a fold was measured and is
  worse than rebuilding per coordinate: 31.04 s and 22537 iterations for
  per-coordinate models against 60.40 s and 25573 for one shared model. The
  infeasible small-mu runs at the end of each coordinate leave a basis the next
  coordinate pays for. One model per coordinate is what the implementation uses.
- Whole-task effect from the measured 78 percent Dantzig share: 3.63x for the
  13.1 factor and 4.08x for the 19.7 factor.
- The advantage grows with `p`, because the frozen path's LP becomes numerically
  harder as the profile Hessian grows. One Dantzig row over the frozen grid on the
  real Hessian of coordinate 1 costs 150.06 s at `p=1000` against 32.46 s for the
  dual form (4.62x) and above 2100 s at `p=2000` against 58.54 s (above 36x). One
  mu-CV call for the same coordinate falls from 684.86 s to 40.12 s at `p=1000`
  (17.07x) and to 22.21 s with the fold split (30.84x), with the same selected
  `mu`, the same candidate vector and defect agreement 1.92e-10. At `p=2000` the
  same call takes 122.77 s serial and 54.55 s over four folds, and the two tables
  agree bit for bit.
- Parity and geometry protection: `tests/testthat/test-dantzig-dual-warm-parity.R`
  (seven blocks, 123 assertions) covers row and mu-CV decisions against the frozen
  path up to `p=200`, the monotone settlement, the failure schema, the dispatcher
  inertness, and the bit-identity of the fold-parallel run against its serial run.
  Measurements and the option search are in
  `COMPUTATIONAL_OPTIMISATION_HIGHS_PATH_REPORT.md` section 9.

### Lever 8: first-stage configuration stack - MEASURED, opt-in
- One stage-1 fit at A02 calls `profile_components_v2` 673 times and takes 90.3 s:
  22 Hessian calls averaging 1.03 s (22.6 s) and 651 score-only calls averaging
  0.104 s (67.6 s). Each call visits all 200 clusters, so the fit performs
  107,600 cluster nuisance solves, which take 51.5 percent of the measured 59.4 s
  of a quieter run of the same fit, at a mean of 0.284 ms and 3.2 criterion
  evaluations per solve. The rest is per-cluster score, objective and bookkeeping
  work.
- Two existing switches act on those halves. `JPDMEQR_CLUSTER_SOLVER=newton`
  replaces the `stats::optim` BFGS nuisance solve with the analytic damped Newton
  path; `JPDMEQR_ASSEMBLY=batch` assembles the Hessian with whole-matrix products.
  One fit, one session: default 60.45 s, `newton` 41.23 s (1.47x), `batch` 50.69 s
  (1.19x), both 40.31 s (1.50x). All four keep `converged=TRUE`,
  `kkt_normalized=9.977e-05` and 32 iterations; `max|dbeta|` is 8.3e-10 for
  `newton` and 1.6e-10 for `batch`. At A08 (`p=2000`) the pair falls from 544.03 s
  to 417.97 s, 1.30x, with the same `kkt_normalized`, the same 28 iterations and
  `max|dbeta| = 1.64e-09`.
- End to end at the same cell, one `PROFILE-DQR` replication with all four
  switches (add `dual_warm` and `cv_cores=4`): 98.0 s and 100.0 s against 569.6 s
  and 651.6 s under the frozen defaults, a whole-task factor of 5.8 to 6.5, with
  the mu-CV part falling from 445.4 s to 31.1 s. Same `mu` on all four target
  coordinates, same `mu` source, same fold counts, relative differences in the
  inference objects below 7e-07. The strict geometry gate passes in the stacked
  configuration (760 cases, 20 seeds).
- The same end-to-end comparison on the two larger estimation cells: at A05
  (`p=1000`, `n=200`) 3825.0 s against 220.8 s, a whole-task factor of 17.32, with
  the same four coordinate rows, the same selected `mu` on each, the same fold
  counts and relative inference differences below 1.1e-08; at A08 (`p=2000`,
  `n=200`) 729.2 s against more than 7.7 h, above 38, where the frozen side is
  composed from its measured stage-1 fit and the stopped mu-CV reference. The factor
  grows with the dimension because the Dantzig share of the frozen runtime does.
- Iteration-loop Hessian refreshes no longer build the four within-fold Hessians
  that only the mu cross-validation reads; the retained fit object still carries
  them. A same-session controlled comparison of one stage-1 fit under `newton` and
  `batch` gives 90.36 s against 73.26 s at A02 (1.233) and 422.44 s against
  356.94 s at A08 (1.184); the 26 A08 Hessian calls take 259.05 s instead of
  322.18 s. Rerun replications at A05 and A08 reproduce the pre-change `coordinates`
  and `theory` tables bit for bit, with only `runtime_sec` differing.
- Configuration provenance: `solver_configuration_v2()` records the solver
  switches, the fold worker count, the first-stage refresh interval and the
  shared-first-stage flag; the
  geometry manifest, the freeze manifest
  and each run directory carry it, and `validate_geometry_gate_v2` accepts a
  manifest only when its record equals the configuration of the current process.
- The first-stage curvature refresh interval is a switch,
  `JPDMEQR_HESSIAN_REFRESH`, default 5. On one stage-1 fit under `newton` and
  `batch` the value 20 gives 1.109 at A02, 1.106 at A05 and 1.279 at A08, with
  every fit accepted inside its `1e-3` residual tolerance and `max|dbeta|` below
  `4e-10`; 40 is slower than 20 at the two smaller cells. The switch reaches the
  two first-stage fits only: the round-5 refit accepts at `1e-7` and keeps the
  frozen value, which is worth 7.1 percent of one A02 replication against 2.4
  percent when the refit is lengthened too. One replication per cell, both sides in
  one session with `cv_cores=2`, gives 156.4 s against 146.0 s at A02, 258.1 s
  against 250.7 s at A05 and 908.4 s against 863.7 s at A08, with the selected
  `mu`, the `mu` source and the fold counts identical at all three cells.
- The lower end of the frozen `mu` grid is inert at the measured cells. A candidate
  reaches the selector only when at least two folds return a finite held-out defect,
  and the `0.02` and `0.05` multipliers never do at A02, A05 or A08: the selector
  returns four or five rows for a grid of eight, and every coordinate selects the
  `0.5` or the `1` multiplier. Removing those two entries leaves the `coordinates`
  table bit-identical at all three cells and leaves the mu-CV term unchanged at A02
  (44.4 s against 45.5 s) and A05 (108.7 s against 108.4 s), while at A08 it falls
  from 409.5 s to 322.4 s and the replication from 829.1 s to 742.8 s, reproduced
  with the runs in the opposite order. Since the gain appears at `p=2000` alone and
  the freeze rule forbids changing a tuning grid without a new versioned registry,
  the eight multipliers stay as they are; adoption would rest on a scan of every
  stored coordinate row showing that the selected `mu` never sat below the new floor.
- The dense share of a switched replication at `p=2000` is large enough to justify
  timing the server BLAS. A stage-1 Hessian assembly is purely R-level dense work,
  and at A08 sixteen of them take 193.4 s of 863.7 s (22 percent), before counting
  the score products, the round-5 refit and the jackknife. This machine's reference
  BLAS moves 2.15 GMAC/s on a 6000x2000 `crossprod` and 1.25 GMAC/s on a 2000-cubed
  `gemm` across 8 logical cores, so a threaded vectorised BLAS plausibly takes
  1.25-1.4x off whole-task runtime at that size. The extrapolation rests on the
  measured share and throughput, and the server should time one replication with and
  without the library. A threaded BLAS makes rounding depend on the thread count, so
  the parity tests and the geometry gate belong in a run with the thread count pinned
  and recorded.
- Parity protection: `tests/testthat/test-solver-configuration.R` covers the
  defaults, the environment-over-option order, the worker count, the refresh
  interval, the manifest comparison and the JSON round trip, and
  `tests/testthat/test-hessian-refresh-switch.R` covers the resolution order, the
  first-stage injection, the isolation of the plain fit and the error on a value
  that is not a positive integer. The first-stage
  switches themselves stay covered by `test-assembly-modes.R`,
  `test-profile-components-identity-opt-in.R` and the stacked gate run. The gate
  script compares the two assembly modes on six deterministic cases (full Hessian,
  four fold Hessians, identity matrix) and records `assembly_modes_checked` and
  `assembly_pass`; `tests/testthat/test-geometry-gate-manifest.R` asserts that the
  gate decision rejects a manifest without that comparison.

## 3. Frozen invariants (must not change)

- Dantzig program min ||omega||_1 s.t. ||H omega - e_k||_inf <= mu;
- mu grid values and the smallest-feasible rule; the 2/4-fold CV structure and
  its loss; all seeds; all registry rows/config/checkpoint versions;
- first-stage formula, bandwidth h_est/h_inf, refit set rule (round 5);
- score/Hessian/sandwich formulas; gate thresholds and V1 semantics.

## 4. Costs

- Implementation, parity suites and geometry re-validation for the adopted
  switches are done locally; one strict gate cycle measures 3.6 minutes on an
  8-core desktop in the stacked configuration (760 cases, 20 seeds), against
  a design estimate of ~1 h under the frozen defaults.
- Every new commit that touches the solver path repeats the gate cycle once,
  in the same configuration as the runs it authorises.
- No effect on the frozen registry contents (same B=200 rows).

## 5. Schedule

1. implemented: Direction 1 and Levers 2, 5, 7, 8 as opt-in switches;
2. parity suites and the strict geometry gate pass on the current tree in the
   stacked configuration;
3. server power-on: fetch assets, install OpenBLAS (Lever 3) and time one
   replication with and without the threaded BLAS;
4. freeze preflight (b200 config) -> launch registry at ~64 cores:
   **~3-4 days wall** (vs 10-12 without the plan);
5. application layer + manuscript backfill in parallel.

Solver substitution and intra-replication fit sharing are implementation
layer; no statistical object changes; the amendment file notes this plan.