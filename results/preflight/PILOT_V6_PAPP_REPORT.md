# PAPP application-matched cells — report (PAPP-HD done; PAPP-LD pending)

## PAPP-HD (n=129, p=500, s=5, tau=0.50, q=1, empirical GSE visit multiset, B=200)

Run: server (Tencent Cloud Ubuntu 26.04, commit 1dcc159 + portable launcher),
36 shards merged; 600 method rows, 2400 coordinate rows, 0 runner errors,
post_refit ok = 1.000, KKT p90 = 6.2e-4.

First stage: l1_lasso med 3.44 -> l1_refit med 2.06; selected_size med 2,
refit_set med 5, omitted_active med 2 (the L1 stage recovers ~3/5 actives at
n=129).

PROFILE-DQR (practical high-dimensional method at the application n):
- coverage 0.330-0.760; SE/SD 0.44-0.60; bias x00002 +0.146.
POSTREFIT-EXACT-H (data-selected exact row, no Dantzig):
- coverage 0.940/0.520/0.935/0.950; SE/SD 0.81-1.14; the under-coverage is
  localised to x00002 (bias +0.167), the coordinate whose active signal is
  missed by selection.
TRUE-SUPPORT: coverage 0.795-0.880, SE/SD 0.66-0.83.

Interpretation: at the application's 129 independent subjects the
high-dimensional practical interval is not calibrated (0.33-0.76), primarily
because the first-stage selection misses ~2/5 actives and the sandwich
under-estimates; the exact-row diagnostic calibrates where the active is in
the data-selected refit set. This is the formal empirical record supporting
round-6 Q3: no confirmatory high-dimensional gene-level Wald/FDR at n=129; the
application keeps the exploratory HD layer + the prespecified low-dimensional
layer (PAPP-LD, pending).

## PAPP-LD (n=129, d=20, s=5, empirical multiset, B=200 per tau in
{0.25, 0.50, 0.75}; ordinary sandwich vs delete-one-cluster jackknife)

Run: cloud server (48 shards, checkpoint-capable runner post-d40d6a9); 2400
rows, 200 replications, status ok = 1.000, no failures.

Coverage (sandwich, mcse 0.009-0.017): tau=0.25: 0.970-0.985;
tau=0.50: 0.935-0.980; tau=0.75: 0.935-0.960. SE/SD 1.02-1.22 (mild
over-dispersion). All 12 (tau x coordinate) cells inside or within mcse of the
[0.88, 0.98] band.

Delete-one-cluster jackknife coverage: tau=0.25: 0.960-0.985;
tau=0.50: 0.935-0.975; tau=0.75: 0.930-0.960; SE/SD 1.00-1.19. Per-cell
sandwich - jackknife coverage difference <= 0.01 in 11/12 cells (max 0.01).

Interpretation: in the prespecified low-dimensional layer (d=20 < n/log n at
n=129) the ordinary unsmoothed cluster sandwich is calibrated, and the full
delete-one-cluster jackknife confirms it (no material disagreement). This is
the empirical basis for the application's two-layer announcement: the
high-dimensional layer at n=129 is exploratory only (PAPP-HD: 0.33-0.76), the
prespecified low-dimensional layer is inferential (PAPP-LD: 0.93-0.99).