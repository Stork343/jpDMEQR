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
   OpenBLAS installed (Lever 3 in COMPUTATIONAL_OPTIMISATION_PLAN.md).
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

Accept: every asset written with `implementation_commit = 320bb43…`,
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
`commit_sha == 320bb43…` and all five booleans true.

## 4. Strict geometry (Stage 2)

```bash
Rscript --vanilla scripts/simulation/00_validate_profile_geometry.R --strict --seed-count 20
```
Accept: `results/preflight/geometry/manifest.json` has `pass=true`,
`strict=true`, `seeds=20`, `commit == 320bb43…`, dantzig checks present.

## 5. Micro-preflight (Stage 5)

```bash
Rscript --vanilla scripts/simulation/02_run_main.R --config=config/simulation_preflight.csv
```
Accept: `results/preflight/micro_preflight_manifest.json` has `pass=true`,
`commit == 320bb43…`, zero runner failures.

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
`commit == 320bb43…`, mechanism evidence present. (The se_ratio/coverage
columns are V1 diagnostic bands, not blocking criteria.)

## 7. Final freeze preflight (Stage 7) — THE GATE

```bash
Rscript --vanilla scripts/simulation/05_preflight_freeze.R --config=config/simulation_main_b200.csv
```
This writes `results/preflight/freeze_manifest.json` with
`config_sha256` of the **B=200** registry (NOT `simulation_main.csv`, which
contains unauthorised B=500/1000 replication counts). It stops with a clear
error unless governance + registry + benchmark + geometry + targets + micro +
pilot all pass against commit 320bb43….

Accept: `overall: TRUE`.

## 8. Final run (Stage 8) — B=200 only

```bash
# jobs = shard-level parallelism, JPDMEQR_CV_CORES = in-task mu-CV parallelism.
# On the paid server (48 logical cores): e.g. jobs=12, cv_cores=4.
JPDMEQR_CV_CORES=4 Rscript --vanilla scripts/simulation/02_run_main.R \
  --config=config/simulation_main_b200.csv --jobs=12
```
`02_run_main.R` refuses to start unless the freeze manifest matches the
current commit AND the config checksum (verify_freeze_manifest_v2). A refusal
is a correct outcome — investigate, never use `--development` for manuscript
results.

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