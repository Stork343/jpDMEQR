# Lever-1 (chained OSQP) adjudication — measured decision incompatibility

Date: 2026-09-03 (local machine, R 4.6.1, osqp 1.0.0, CVXR/CLARABEL)
Status: **Lever 1 must NOT be enabled as written.** Decision-incompatible with
the frozen smallest-feasible-mu rule; no implementation-layer parameter choice
fixes it. Report to study lead for re-planning.

## 1. Context

`docs/COMPUTATIONAL_OPTIMISATION_PLAN.md` (2026-09-01) proposed replacing the
CVXR/CLARABEL dense LP solves in the Dantzig mu-CV (4 coords x 8 mu x 4 folds =
128 LPs per task, ~75-80% of task time) with a persistent OSQP model and
warm-start chaining, claiming 20-50x on the dominant term ("first feasible mu
in 0.8s including setup"). The implementation exists (commit e1263d8) and is
selected by `JPDMEQR_DANTZIG_SOLVER=osqp`, but defaults to CVXR/CLARABEL.

## 2. Measured facts (this machine, correlated design H)

Single-row solve (`solve_dantzig_row_v2`, mu grid 10^-1.6..10^0.5):

| p   | CLARABEL | OSQP (current settings) | notes |
|-----|----------|-------------------------|-------|
| 200 | 1.27 s   | 9.14 s                  | CLARABEL returns after 1 attempt (first feasible mu) |
| 500 | 6.75 s   | 120.9 s                 | OSQP solves all 8 mu, no early exit |

mu-CV selection path (`select_mu_inverse_hessian_cv_v2`, 4 folds x 8 mu):

| p   | CLARABEL | OSQP repaired (eps=1e-6, polish) | speedup | selected mu |
|-----|----------|----------------------------------|---------|-------------|
| 200 | 10.6 s   | 9.0 s                            | 1.2x    | 0.0251 vs **0.2** |
| 500 | 126.0 s  | 66.8 s                           | 1.9x    | 0.0251 vs **0.398** |

Parameter scan on the smallest mu (p=200, mu=0.02512, the mu where CLARABEL is
first feasible):

| settings | result |
|----------|--------|
| CLARABEL | residual 0.025119 <= mu + 1e-6 -> feasible (1 LP, 1.3 s) |
| OSQP eps=1e-8 polish 20k iters | residual 0.025155 -> INFEASIBLE |
| OSQP eps=1e-8 polish 100k iters | residual 0.025163 -> INFEASIBLE (14.8 s) |

OSQP first feasible mu on the grid at best settings: **0.0501** vs CLARABEL
**0.0251**.

## 3. Root cause

The frozen Dantzig row is an L1-constrained LP whose feasibility boundary
(`||H omega - e_k||_inf <= mu`) is *exactly tight* at the smallest feasible mu
(residual ~ mu). OSQP is an ADMM first-order method: on such tightly binding
L1 problems it stalls at residual ~ mu + 4e-5 and does not converge within
100,000 iterations, while CLARABEL (interior point) resolves the same problem
in one solve. This is a method-level limitation, not a tuning issue.

Additionally the current chain implementation solves all 8 mu without early
exit (CLARABEL stops at the first feasible mu), compounding the slowdown;
fixing that alone still leaves decision incompatibility at the tight mu.

## 4. Consequences

- Enabling `JPDMEQR_DANTZIG_SOLVER=osqp` would change the frozen mu selection
  (0.0501 vs 0.0251) and therefore the debiased estimates / intervals;
  prohibited by the freeze and by AGENTS.md (no decision change).
- The plan's 20-50x claim is not reproduced (1.2-1.9x measured on the CV path
  after repair). The plan's "0.8 s" measurement was not decision-parity
  validated and does not match this machine's 9-120 s.

## 5. Recommended replacement directions (implementation layer only)

1. **Parallelise the 32 independent CLARABEL solves inside the mu-CV** (each
   (fold, mu) solve is independent; aggregate defects in fold order -> results
   bit-identical, decision unchanged). Wall-clock for the dominant term
   divided by available cores per task; combine with coarse task-level
   parallelism to avoid nested-fork memory blowup.
2. **OpenBLAS / threaded BLAS on the server** (2-4x on all dense ops,
   including CLARABEL's linear algebra), as planned.
3. **Lever 2: share the calibrated first stage across methods within a
   replication** (PROFILE-DQR / POP-H / TRUE-SUPPORT currently refit the same
   calibrated first stage; bit-identical, 1.2-1.5x on multi-method cells).
4. **Shard-level load balancing**: replication cost varies by ~100x across
   experiments (p=200 vs p=2000); the wall clock of a shard is driven by its
   slowest task — finer-grained dynamic scheduling reduces tail latency.
5. Re-measure CLARABEL's fixed cost: CVXR problem construction per solve may
   be avoidable (call CLARABEL's native interface), worth profiling.

## 5a. Direction 1 implemented (2026-09): parallel mu-CV CLARABEL grid

Selection accepted by the study lead. Implemented in `select_mu_inverse_hessian_cv_v2`
(CLARABEL branch; `cv_cores <= 1` never leaves the frozen serial path) with
top-level worker `dantzig_cv_solve_task_v2` and dispatcher
`dantzig_cv_clarabel_parallel_v2` in `R/v2_profile_v2_part04.R`.

- Bit-identity: the (fold x mu) grid is solved concurrently and aggregated in
  the original (mu outer, fold inner) order; each solve is a pure function of
  its inputs, so mean/sd, the argmin-tie rule, and the infeasible-fold drop
  rule see exactly the serial numbers.
- Backends: Linux fork-based `mclapply` (keep `jobs` at 1, use shards —
  no nested fork); Windows deterministic PSOCK with `parLapplyLB` dynamic load
  balancing (CLARABEL per-cell wall time varies 2.9-7.4 s at p=500; static
  round-robin would strand each fold's whole mu chain on one worker).
- Control: option/env `jpDMEQR.cv_cores` / `JPDMEQR_CV_CORES`, default 1
  (serial, bit-identical); wired into `launch_sharded_pilot.R` and
  `02_run_main.R`.
- Measured (8-core Windows desktop, p=500, 4 workers): serial 136.8 s vs
  parallel 49.7 s per coordinate-CV = 2.75x; selected mu, min defect, and
  every candidate row bit-identical. p=120: 1.93x.
- Gate: new `tests/testthat/test-cv-parallel-parity.R` asserts bit-identity
  serial vs cv_cores=2/4 across coordinates; full v2 test suite green.
- Expected registry effect at fixed core count: dominant-term wall time
  divides by cv_cores (e.g. jobs=12 x cv_cores=4 on 48 logical cores) with
  no statistical-object change.

## 6. Verification evidence

- `tmp_dantzig_bench.csv`, `tmp_cv_bench.csv`, `tmp_osqp_repaired_cv.csv`
  (benchmark outputs; kept outside the repo).
- Reproduction scripts `tmp_osqp_scan.R`, `tmp_osqp_final.R`;
  Direction-1 scripts `tmp_cv_parity.R`, `tmp_cv_bench_p500.R`,
  `tmp_cv_breakdown.R`, `tmp_cv_lb_compare.R`, `tmp_psock_probe.R`
  (all outside repo; remove after study-lead review if desired).
- Existing parity unit test `tests/testthat/test-dantzig-solver-parity.R`
  covers p=20/60/150 only and did not catch the p>=200 CV decision mismatch —
  extend it before any solver substitution is ever re-attempted.