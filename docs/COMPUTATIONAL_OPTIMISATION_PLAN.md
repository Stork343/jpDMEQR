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

### Lever 1: persistent OSQP model with warm-start chaining (the big one)
- The 128 CV solves per task share structure: per fold, the constraint matrix A
  and objective P are IDENTICAL across all mu and all 4 coordinates; only the
  box bounds (l,u) change (e_k and mu enter the bounds only).
- Realisation: osqp's C-level API (osqpSetup/osqpUpdate/osqpWarmStart,
  exposed inside the installed osqp 1.0.0 namespace) keeps ONE KKT
  factorisation per fold; every subsequent solve costs O((n+m)^2) with warm
  starts.
- Measured locally: p=500, first feasible mu in 0.8 s including setup
  (CLARABEL: 3.5-4.3 s per single solve); the 4-fold CV is projected at
  10-20 s/task vs ~373 s -> 20-50x on the dominant term. At p=2000 the
  one-per-fold factorisation is amortised over 32 solves/coord-grid ->
  ~30x on those cells.
- Parity: same feasibility acceptance (residual <= mu + 1e-6 and solver
  status), same first-feasible-mu rule; decision parity vs CLARABEL verified
  by the existing solver-parity test and re-checked by the strict geometry
  Dantzig validation in the freeze pipeline (the pipeline runs ONE solver
  consistently).
- Expected registry effect: ~29,500 core-h -> ~6,000-8,000 core-h.

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
  (first stage Hessian assembly, OSQP factorisation steps).
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