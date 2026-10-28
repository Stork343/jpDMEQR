# Server readiness runbook (pre-launch checklist)

**Status:** DRAFT — produced by the onshore readiness audit (2026-09-03) at commit
`320bb43c0accaa99b6a48f2099a99cbb622bfe1f` on branch `gpt/simulation-freeze-resolution`.

**Purpose:** the exact, ordered sequence that must be run on the paid server
BEFORE the final B=200 registry is started. Each step carries its own gate;
step 8 (freeze preflight) is the single switch that authorises the main run.
If any step fails, stop and report — do not bypass.

---

## 0. Preconditions (outside the repository)

1. **GitHub reachable from this machine.** The local branch is **4 commits
   ahead** of `origin/gpt/simulation-freeze-resolution` (0183d3b...320bb43).
   Push first:
   ```bash
   git push origin gpt/simulation-freeze-resolution
   ```
   As of 2026-09-03 the proxy to github.com:443 was failing; nothing can be
   pulled or pushed until that is fixed. This is the only external blocker.
2. **Server software:** R 4.6.1 (same as development machine), packages:
   `CVXR, CLARABEL, SCS, OSQP, quantreg, digest, jsonlite, testthat, pkgload`.
   **OpenBLAS: measure it on the server, and expect more than the earlier
   estimate.** The dominant term under the frozen defaults is CLARABEL, a Rust
   solver using its own `qdldl` factorisation through a self-contained DLL, so it
   does not pass through R's BLAS. The switched stack changed the composition:
   with `dual_warm` the mu-CV term is HiGHS in a DLL, while at A08 one replication
   spends 193.4 s of 863.7 s (22 percent) in the first-stage Hessian assembly
   alone, all of it in R-level dense products, with the score products, the
   round-5 refit and the jackknife adding further dense work. The development
   machine's reference BLAS moves 2.15 GMAC/s on a 6000x2000 `crossprod` (11.17 s)
   and 1.25 GMAC/s on a 2000-cubed `gemm` (6.42 s) on an 8-core box, a few percent
   of what the cores could reach with vectorised threaded kernels. A threaded
   OpenBLAS therefore plausibly takes **1.25-1.4x off whole-task runtime at
   `p=2000`**, concentrated in the first stage. That figure extrapolates from the
   measured share and throughput, so time one replication with and without it
   before scheduling around it. A different BLAS reorders summations and moves the
   last bits of every product, and a threaded BLAS makes a product's rounding
   depend on the thread count, so run the serial-versus-parallel parity tests and
   the geometry gate with the thread count pinned (for example
   `OMP_NUM_THREADS=1`) and record the library, its version and that setting in the
   session/hardware manifest.
3. **Repository at the exact commit:** `320bb43c0accaa99b6a48f2099a99cbb622bfe1f`
   (do not run on any other commit — every gate below compares against it).
4. Working directory = repository root; do not place under a synced folder.

## 1. Source integrity (Stage 0)

```bash
Rscript --vanilla -e "for (f in list.files('R', pattern='^v2_.*\\.R$', full.names=TRUE)) invisible(parse(f)); cat('parse ok\n')"
Rscript --vanilla -e "if (!requireNamespace('pkgload', quietly=TRUE)) stop('pkgload missing'); pkgload::load_all('.'); testthat::test_dir('tests/testthat', reporter='summary')"
```
Accept: all v2-related tests green; the 7 known legacy failures
(generate_mixed_qr_data / Simulation_Study / default_lambda_clime) are
pre-existing and unrelated to the v2 pipeline.

## 2. Target assets (Stage 4) — all assets, in one run

`results/intermediate/` is gitignored, so the server has **none** of the
profile / audit / population-direction assets. Rebuild everything on the
current commit in one call (this is what also fixes the stale
`implementation_commit` on the 23 POP-H assets):

```bash
Rscript --vanilla scripts/simulation/04_build_profile_targets.R \
  --config=config/simulation_main_b200.csv
```
Expect: 5 profile targets (C06,C07,C14,C15,C16), 4 structural audits
(E03,E04,E05,E06), 23 population-direction assets (A01–A12, C06, C07, C14,
C15, C16, D01–D06), each with 100000 clusters × 4 repeats (final scale).
Runtime on 48 logical cores: roughly 1–2 h. On Windows set
`JPDMEQR_REPEAT_CORES`/PSOCK as documented in 04; on Linux the script uses
mclapply automatically.

Accept: every asset written with `implementation_commit = the current commit`,
`max_nuisance_gradient <= 1e-7`, direction residuals <= 1e-5, and a fresh
`results/intermediate/.../build_manifest.json`.

## 3. Benchmark acceptance (Stage 3)

```bash
Rscript --vanilla scripts/benchmarks/accept_patch_adapters.R       # PROFILE-DQR, POOLED-QR-LASSO, LQMM, PROFILE-DQR-TRUE-SUPPORT
Rscript --vanilla scripts/benchmarks/accept_remaining_benchmarks.R  # BIAS-ADJ-LQMM, DOUBLE-PEN-QLMM, QGEE-SCAD, QIF-SEE
Rscript --vanilla scripts/benchmarks/accept_oracle_adapters.R       # PROFILE-DQR-TRUE-NUISANCE, PROFILE-DQR-POP-H, PROFILE-DQR-SPLIT
Rscript --vanilla scripts/benchmarks/accept_sqr_debiased_iid.R      # SQR-DEBIASED-IID
```
Accept: all 12 `results/preflight/benchmarks/<METHOD>/acceptance.json` have
`commit_sha == the current commit` and all five booleans true.

## 4. Strict geometry (Stage 2)

```bash
Rscript --vanilla scripts/simulation/00_validate_profile_geometry.R --strict --seed-count 20
```
Accept: `results/preflight/geometry/manifest.json` has `pass=true`,
`strict=true`, `seeds=20`, `commit == the current commit`, dantzig checks present
(`dantzig_checked`, `dantzig_pass`), and the assembly-mode comparison present
(`assembly_modes_checked`, `assembly_pass`; six deterministic cases comparing the
full Hessian, the four fold Hessians and the identity matrix across the two
assembly modes).

## 5. Micro-preflight (Stage 5)

```bash
Rscript --vanilla scripts/simulation/02_run_main.R --config=config/simulation_preflight.csv
```
This is the one registry that runs without a freeze manifest, because its own
manifest is an input to that freeze (`requires_freeze_manifest_v2`). The
`profile_mc` cells of this registry need their target assets built on the same
commit as the run (`PF04`, `PF05`, `PF06`), and `PF04` also needs its
population-direction asset.
Accept: `results/preflight/micro_preflight_manifest.json` has `pass=true`,
`commit == the current commit`, zero runner failures.

## 6. Pilot (Stage 6) + pilot gate

`pilot_gate.json` must be re-evaluated on the current commit. The pilot raw
runs are gitignored, so either (a) re-run the frozen pilot on the server, or
(b) copy the local `results/raw/simulation/pilot_v5_B200*` across and
re-evaluate:

```bash
# (a) fresh pilot on the server:
Rscript --vanilla scripts/simulation/01_run_pilot.R --config=config/simulation_pilot_v5.csv --run-id=pilot_v5_B200 --max-reps=200 --jobs=<N>
# then evaluate:
Rscript --vanilla scripts/simulation/03_summarise.R --run-id=pilot_v5_B200 --pilot-gate=true
```
Accept: `results/preflight/pilot_gate.json` has `pass=true`,
`commit == the current commit`, mechanism evidence present. (The se_ratio/coverage
columns are V1 diagnostic bands, not blocking criteria.)

## 7. Final freeze preflight (Stage 7) — THE GATE

```bash
Rscript --vanilla scripts/simulation/05_preflight_freeze.R --config=config/simulation_main_b200.csv
```
This writes `results/preflight/freeze_manifest.json` with
`config_sha256` of the **B=200** registry (NOT `simulation_main.csv`, which
contains unauthorised B=500/1000 replication counts). It stops with a clear
error unless governance + registry + benchmark + geometry + targets + micro +
pilot all pass against the current commit.

Accept: `overall: TRUE`.

## 8. Final run (Stage 8) — B=200 only

```bash
# jobs = shard-level parallelism, JPDMEQR_CV_CORES = in-task mu-CV parallelism,
# JPDMEQR_DANTZIG_SOLVER = which path solves the frozen Dantzig LP (see below).
# On the paid server (48 logical cores): e.g. jobs=12, cv_cores=4.
JPDMEQR_DANTZIG_SOLVER=clarabel_native JPDMEQR_CV_CORES=4 \
  Rscript --vanilla scripts/simulation/02_run_main.R \
  --config=config/simulation_main_b200.csv --jobs=12
```
`02_run_main.R` refuses to start unless the freeze manifest matches the
current commit AND the config checksum (verify_freeze_manifest_v2). A refusal
is a correct outcome — investigate, never use `--development` for manuscript
results.

### Performance knobs (wall-clock / cost only — never statistical objects)

| knob | effect | notes |
|---|---|---|
| `jobs` | shard / task-level parallelism | multiplies throughput; does not reduce total core-hours |
| `JPDMEQR_CV_CORES` | parallelises the mu-CV Dantzig grid (the 75-80% term) | default 1. With `clarabel_native` the (fold, mu) cells go to workers: measured **4.2-5.3x** on that term at 4 cores. With `dual_warm` the unit is a whole fold (4 in the frozen rule): measured **3.42x and 4.02x** against its own serial run at 4 cores; memory, not cores, is the limit at `p=2000` (see the note under the table) |
| `JPDMEQR_DANTZIG_SOLVER` | which solver path solves the frozen LP | default (unset) = frozen CVXR path. `clarabel_native` = compact native CLARABEL form, measured **1.30-1.64x** on that term. `dual_warm` = dual form on a persistent HiGHS simplex model with the basis carried along the mu grid, measured **4.05x and 6.00x** serial on the same two coordinates of A02, **13.1x and 19.7x** with the fold split at 4 cores; the factor grows with `p`, reaching **17.07x serial and 30.84x over four folds** at A05 (`p=1000`) and above **55x and 125x** at A08 (`p=2000`), where one frozen candidate alone costs 428.08 s against 3.17 s. `highs` = the compact form solved by HiGHS, **not recommended** (adjudication below) |
| `JPDMEQR_CLUSTER_CORES` | cluster-level parallelism inside a task | default follows `jobs`. On Linux keep `jobs=1` to avoid nested-fork memory blowup |
| `JPDMEQR_CLUSTER_SOLVER` | how the per-cluster nuisance problem is solved in the first stage | default `bfgs` (frozen reference). `newton` uses the analytic gradient and curvature. On one stage-1 fit at A02 the whole stage falls from 60.45 s to 41.23 s, a factor of **1.47**, with `max|dbeta| = 8.3e-10`; paired with `batch` at A08 (`p=2000`) it falls from 544.03 s to 417.97 s, a factor of **1.30**, with `max|dbeta| = 1.64e-09`, the same `kkt_normalized` and the same 28 iterations |
| `JPDMEQR_ASSEMBLY` | how the first-stage Hessian is assembled | default `sequential` (frozen reference, per-cluster accumulation). `batch` keeps the cluster rows and assembles with whole-matrix products. On the same fit the stage falls from 60.45 s to 50.69 s, a factor of **1.19**, with `max|dbeta| = 1.6e-10`; paired with `newton` it gives the 1.50 factor at A02 and the 1.30 factor at A08 quoted in the row above |
| `JPDMEQR_HESSIAN_REFRESH` | how many iterations the first stage goes between curvature refreshes | default 5 (frozen reference). Each refresh costs a `p x p` assembly, so a longer interval buys time at the price of a staler momentum step. `20` measures **1.109** at A02, **1.106** at A05 and **1.279** at A08 on one stage-1 fit, with the fit still converging inside its acceptance level (`kkt_normalized` 2.34e-04, 5.81e-04 and 1.13e-04 against the 1e-03 tolerance) and `max|dbeta|` below 4e-10; `40` is worse than `20` at A02 and A05 because the extra iterations outweigh the saved assemblies. The round-5 refit, the jackknife, the population target assets and the geometry audit keep the frozen 5 |
| `JPDMEQR_SHARED_FIRST_STAGE` | shares the calibrated first-stage fit across methods inside one replication | default off (frozen reference). When on, PROFILE-DQR, PROFILE-DQR-POP-H and POSTREFIT-EXACT-H reuse the same deterministic fit object; outputs are bit-identical to independent fits (`test-shared-first-stage.R`). Saves one stage-1 fit per additional sharing method per replication: 40.31 s at A02 and 417.97 s at A08 under the switch stack. Bound into the geometry manifest like every other knob. See Lever 2 in `COMPUTATIONAL_OPTIMISATION_PLAN.md` |

**Recommended server setting:** the four switches together,
`JPDMEQR_DANTZIG_SOLVER=dual_warm`, `JPDMEQR_CV_CORES=4`,
`JPDMEQR_CLUSTER_SOLVER=newton` and `JPDMEQR_ASSEMBLY=batch`, extended by
`JPDMEQR_HESSIAN_REFRESH=20` (quantified two paragraphs below) and
`JPDMEQR_SHARED_FIRST_STAGE=1`, which drops the duplicated calibrated
first-stage fits in multi-method cells at bit-identical results. On the registered
cell A02 one `PROFILE-DQR` replication takes 98.0 s and 100.0 s in two sessions
against 569.6 s and 651.6 s under the frozen defaults, a **factor of 5.8 to 6.5
on whole-task runtime**, with the mu-CV part falling from 445.4 s to 31.1 s and
the first stage from 60.45 s to 40.31 s. The two stages are independent, so that
range holds under either session's machine load. Both settings select the same
`mu` on all four target coordinates, the same `mu` source and the same fold
counts; the inference objects agree to 7e-07 in relative terms. The fold split is
bit-identical to the serial run of the same call.

Two of the four switches act on the first stage, and the code now also drops the
within-fold Hessians that the fit's iteration refreshes never hand to a consumer.
That omission is worth a further **1.233** on one stage-1 fit at A02 and **1.184** at
A08 (`p=2000`), measured by running the same fit twice in one session with the flag
forced on and off; the result tables are unchanged at both sizes. Adding
`JPDMEQR_HESSIAN_REFRESH=20` to the set is worth **1.071** at A02, **1.030** at A05
and **1.052** at A08 on a whole replication, measured the same way, with the
selected `mu`, the `mu` source and the fold counts identical at the three cells.

The gain grows with the dimension, which is where the grid spends its time. At A05
(`p=1000`, `n=200`) one replication takes **3825.0 s under the frozen defaults and
220.8 s under the switches, a factor of 17.32**, with the same four coordinate rows,
the same selected `mu` on each of them, the same fold counts and relative
differences below 1.1e-08 in the inference columns. At A08 (`p=2000`, `n=200`) the
switches take 729.2 s and converge with 28 iterations; the frozen path at that size
needs more than 7.7 h, composed from its measured stage-1 fit of 544.03 s and more
than 6829 s for a single mu-CV coordinate, so no run uses the frozen path at
`p=2000`. With `clarabel_native` in place of `dual_warm`, pair it with
`JPDMEQR_CV_CORES=4` for ~4-5x on the dominant term and roughly **2.7x on whole-task
runtime**.

The configuration a gate was measured under travels with its manifest from here
on: the geometry manifest, the freeze manifest and each run directory record the
solver switches, the worker count, the first-stage refresh interval and the
shared-first-stage flag, and the
freeze preflight accepts a geometry manifest only when that record equals the
configuration of the running process. Adopting any of these switches therefore
requires a geometry gate run in the same configuration.

At `p=2000` the CV worker count runs into memory before it runs out of cores. Every
fold worker builds its own dual model, whose constraint matrix has `2p` rows of `2p`
columns (128 MB at `p=2000`), holds its own basis and, on Windows, receives its fold
Hessian through PSOCK serialisation. On an 8-core, 31.5 GB machine one A08
replication with `cv_cores=4` ended in `precision_solver` with `std::bad_alloc`
after 477.7 s while the same configuration completed in an earlier session. Use
`cv_cores=2` at `p=2000` on a 32 GB box, or confirm free memory before launch;
memory per worker grows as `p^2`.

### The dual form on a persistent simplex model

`JPDMEQR_DANTZIG_SOLVER=dual_warm` solves the dual of the frozen row program,
`min e_k'w + mu||w||_1` subject to `||H'w||_inf <= 1`, on one persistent HiGHS
simplex model per fold and coordinate. The dual's constraint matrix depends on
`H` alone, so the eight candidates of the grid are reoptimisations of a single
model with the basis carried between them; a candidate that returns an unbounded
certificate settles every smaller candidate without a solve. On A02 one mu-CV
call takes 28.53 s (coordinate 1) and 18.92 s (coordinate 3) against 115.55 s and
113.49 s under the frozen reference; with `JPDMEQR_CV_CORES=4` the same calls take
8.81 s and 5.75 s. All four runs select the same `mu` (`4.3188077227e-01`), and
the parallel tables equal the serial tables bit for bit. Full measurements, the
verification test file and the option search are in
`COMPUTATIONAL_OPTIMISATION_HIGHS_PATH_REPORT.md` section 9.

### The HiGHS path is adjudicated inadvisable

HiGHS wins on a well-conditioned construction of the same LP and loses on the
real profile Hessian. The comparison at `p=500`, `n=200`, the registered cell
A02, one Dantzig row over the frozen grid:

| Hessian | condition | CLARABEL compact | HiGHS |
|---|---|---|---|
| synthetic construction | 6.04 | 4.9 s / 5.5 s | 3.3 s / 2.6 s |
| real profile Hessian | 1.604e+19 | 15.8 s / 15.0 s | 181.3 s / 181.6 s |

Both paths select the same `mu` (1.079702e-01) on the real Hessian, so the
decision is unchanged while the cost inverts. The real profile Hessian at this
cell is numerically singular: smallest singular value 2.589e-19 against a
largest of 4.153e+00. The compact CLARABEL path itself reports
`almost_primal_infeasible` and one `numerical_error` on the candidates before
the selected one; the HiGHS run hit its 60 s per-solve limit on the first three
of four candidates. Adjudication: `highs` stays in the code behind its switch as
instrumented evidence, and no server run uses it.

### Provenance requirement when `clarabel_native` is enabled

The compact form reproduces the frozen CVXR reference at machine precision
(max `|d omega|` ~2.2e-16) and, decisively, gives **identical selected mu and
identical fold-mean defect (0.00e+00 difference)** at p=200 — verified end-to-end
in `tmp-opt/verify_cv_endtoend.R` and covered by
`tests/testthat/test-dantzig-native-parity.R` at p = 20/60/150/200. It is
therefore decision-identical, but **not byte-identical**, to the reference.

On the registered cell A02 (`n=200`, `p=500`, the numerically singular profile
Hessian) the factor is smaller than on the well-conditioned constructions: one
mu-CV call takes 99.6 s under the frozen reference and 85.4 s and 67.3 s under
the compact form in two separate runs, so the factor on this cell lies between
1.17 and 1.48 and moves with machine load. A replication issues four such calls,
so the saving is roughly 60 to 130 s per replication at that cell. The 1.30-1.64x
in the knob table comes from the earlier measurement set and is not what this
cell shows.

Consequences for the record:

- the freeze manifest's commit pin must be advanced to the commit that contains
  this change (the runbook's `320bb43…` pins are stale once these files are
  committed);
- the run's session/provenance record must state which Dantzig path was used, so
  that the numeric output can be attributed;
- if a byte-identical reproduction of the `320bb43…` reference is needed, simply
  leave `JPDMEQR_DANTZIG_SOLVER` unset — the default path is unchanged.

---

## What the onshore audit already confirmed (do not redo)

- 29 v2 R files parse clean on the current commit (Stage 0).
- `validate_registry_contract_v2(root, config, "main")` is used by the freeze
  preflight for **any** config, so the B=200 registry is validated under the
  main-registry contract.
- The 5 profile_target + 4 structural_audit asset contents have matching
  dependency hashes; the only reason they are rebuilt is that the server
  cannot see the gitignored local copies — rebuilding on the current commit
  satisfies the final gate trivially.
- The old `freeze_manifest.json` fails because it was generated on commit
  93d1362 against `simulation_main.csv` (78 rows incl. unauthorised B=500/1000
  cells); it is superseded by step 7.
- 15 `results/preflight/*` JSON files carry local edits (evidence re-stamped
  to b9534f3) that are **not** to be committed — the server run overwrites
  them with 320bb43-stamped versions, which are then committed by the runner.

## Failure handling

- Any step ≠ pass: stop, save the log, report back. Do not edit the frozen
  registry, seeds, methods, grids, or thresholds (AGENTS.md §8.3).
- A B=200 cell that fails to converge is a reportable outcome, not a reason
  to change the code mid-run.
- Performance knobs (jobs, cv_cores) change only wall-clock time, never the
  frozen statistical objects.

## GitHub connectivity: environment-proxy gotcha (observed 2026-09-03)

Symptom on the development machine:

```text
fatal: unable to access 'https://github.com/Stork343/jpDMEQR.git/':
Failed to connect to github.com:443 over proxy 127.0.0.1 after ~2 s:
Could not connect to server
```

Diagnosis: the machine's `HTTP_PROXY` / `HTTPS_PROXY` (and lowercase)
environment variables point at `http://127.0.0.1:7897`, but **that proxy port
does not actually serve** — while a direct TCP connection to
`github.com:443` succeeds. Git honours the env proxy and fails; setting
`git config http.proxy` or `git -c http.proxy=...` does not fix it because
the env vars take precedence.

Fix (clear the env proxies in the shell where git runs, then retry):

```powershell
$env:HTTP_PROXY=''; $env:HTTPS_PROXY=''; $env:http_proxy=''; $env:https_proxy=''
git fetch origin
git push origin gpt/simulation-freeze-resolution
```

Verification that the proxy (not the network) is the problem:

```powershell
Test-NetConnection github.com -Port 443   # True => network OK, proxy is the culprit
curl.exe -x http://127.0.0.1:7897 -sI --connect-timeout 5 https://github.com  # fails => proxy dead
git -c http.proxy= -c https.proxy= ls-remote origin HEAD   # works => direct route fine
```

This is a machine-environment issue, not a repository issue; if a server
shell shows the same symptom, apply the same env-clear before any
`git pull`/`git push`. After a successful push the branch
`gpt/simulation-freeze-resolution` tracks the exact commit
`745befd…` (runbook + 320bb43 parallelism + prior commits) on the remote.