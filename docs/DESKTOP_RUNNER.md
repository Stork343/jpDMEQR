# Desktop runner kit (offline compute carrier)

**Purpose:** run the heavy R simulation cells on an idle, internet-disconnected
Windows desktop (i7-14700, 32GB RAM, Win10) that has no GitHub access.

## 1. What to copy to the desktop (via USB or local SMB share)

1. **R installation + packages**: from the development machine, copy the whole
   R tree (`E:\Dev\Apps\R`) plus the user library (package lib) to e.g.
   `C:\Dev\Apps\R` on the desktop. R must stay the same version (4.6.1) and
   the library must contain: CVXR, CLARABEL, SCS, OSQP, quantreg, digest,
   jsonlite, testthat, pkgload, parallel (base).
2. **Repository**: either
   - `git bundle create jpdmeqr.bundle gpt/simulation-freeze-resolution`
     (offline-safe), or
   - a plain copy of the working tree at the pinned commit.
   Place it at `C:\jpDMEQR` (NOT under OneDrive/OneDrive-synced folders).
3. **This document**.

## 2. Desktop-side verification (before any long run)

```powershell
cd C:\jpDMEQR
Rscript -e "for (f in list.files('R', pattern='^v2_.*_part[0-9]+\\.R$', full.names=TRUE)) invisible(parse(f)); cat('parse ok\n')"
Rscript -e "pkgload::load_all('.'); testthat::test_dir('tests/testthat', reporter='summary')"
Rscript scripts/simulation/06_variance_ladder_v6.R smoke_ladder P06 ... # or a 1-rep registry smoke:
```

A 1-replicate registry smoke (P05, B=1) must finish with 0 errors and KKT
acceptance before any full cell is started:

```powershell
Rscript <launcher> pilot_desktop_smoke P05 1 1
```

## 3. Launching cells

Windows: one `Rscript` process per shard; size shards by memory
(~2-4 GB/process at p=500-500 turn n=400-800):
recommended **8 shards** on this desktop (32 GB).

Example (PAPP-HD, B=200, 8 shards -> 25 reps each, run from separate Prompts):

```powershell
# shard s in 1..8
Rscript scripts/simulation/01_run_pilot.R --config=config/simulation_papp.csv --run-id=papp_hd_B200 --max-reps=200 --experiments=PAPP-HD --jobs=1 --development
# NOTE: use the sharded launcher (launch_v3_pilot.R) with n_shards=8 for
# replicate-subset sharding; see the launcher's args.
```

Actual command (launcher pattern used by the development machine):

```powershell
Rscript C:\jpDMEQR\scripts\...\launch_v3_pilot.R papp_HD_B200 PAPP-HD 200 1 8 <s>
```

where `<s>` is the shard index 1..8. Run each shard in its own PowerShell
window / background job. Each shard writes
`results/raw/simulation/papp_HD_B200_s<s>/`.

## 4. Copy results back

Copy the whole `results/raw/simulation/<run_id>_s*` directories back to the
development machine (USB/SMB). Merging, analysis, reporting and git commits
happen on the development machine.

## 5. Determinism and provenance

- The repository must be at the pinned commit; the runner records
  `implementation_commit.txt` and `hardware.json` (CPU/RAM/OS) per run, so
  the freeze manifest will truthfully record the desktop hardware.
- R 4.6.1 on both machines + tolerance-based acceptance (KKT <= 1e-3) keeps
  cross-machine numeric differences inside the audit tolerances; the unit
  tests + 1-rep smoke are the gate before any long run.
- No git operations happen on the desktop.

## 6. Target workload for the desktop

1. PAPP-HD (B=200, n=129, p=500, empirical visit multiset) — ~3 methods.
2. PAPP-LD25/50/75 (B=200 each, n=129, d=20, direct unpenalised profile;
   ordinary sandwich vs delete-one-cluster jackknife) — dedicated runner.
3. Final B=500/1000 main registry (only after the freeze manifest passes).
4. Any re-runs that the theory round requires.

Keep OneDrive OFF on the desktop; do not place the repo/results under a
synced folder (file locks and IO contention measurably slowed runs on the
development machine).