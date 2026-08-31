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

Status: running on the cloud server (~456 CPU-h total; first reps complete,
ETA ~21:00). Report section filled after completion.