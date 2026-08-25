# Round-5 paired mechanism probe (P05-v5 / P06-v5, B=50) — report

Runs: results/raw/simulation/pilot_v5_P05_B50 (n=200) and pilot_v5_P06_B50
(n=400), 3 shards each, 0 runner errors, nested clusters (P06 uses P05 seeds;
seed_experiment=P05 in run_request.json). Commit 0dcc159.
Spec: METHOD_SPECIFICATION_ROUND5_AMENDMENT.md. Diagnostic only; no formal
gate from B=50 (PILOT_GATE_THEORY_ACTIONS.txt items 14-15).

## 1. Post-refit mechanism: confirmed

| metric (PROFILE-DQR) | P05-v5 n=200 | P06-v5 n=400 |
|---|---|---|
| l1_error_lasso med | 3.243 | 2.273 |
| l1_error_refit med | 1.569 | 0.246 |
| l2_error_refit med | 0.874 | 0.110 |
| paired l1_refit diff (nested) | — | -1.276 |
| refit_set_size med | 6 | 7 |
| omitted_active_count med | 1 | 0 |
| post_refit ok | 1.000 | 1.000 |
| post_refit gradient med / hess cond med | 6.1e-9 / 6 | - / - |
| KKT normalized med | 1.2e-4 | 1.7e-4 |
| status ok | 1.000 | 1.000 |

The zero-L1 profile refit at h_inf removes the L1 shrinkage bias: at n=400 the
refit error is 0.246 (parametric scale, l2 0.11) with zero omitted actives,
100% refit acceptance, and no dimension-gate failures.

## 2. Four-layer comparison (round-5 diagnostic question answered)

n=200 (PROFILE-DQR / POSTREFIT-EXACT-H / POP-H / TRUE-SUPPORT):
- active coverage: 0.82/0.36 | 0.08/0.00 | 0.96/0.88 | 0.92/0.86
- A_k active: 0.46/1.22 | 0.75/1.30 | 0.63/0.87 | -
- se/sd active: 0.68/0.27 | 1.22/0.49 | 1.03/0.86 | 0.92/0.84

n=400:
- active coverage: 0.86/0.76 | - | 0.98/0.94 | 0.92/0.90
- A_k active: 0.21/0.27 | - | 0.12/0.15 | -
- se/sd active: 0.75/0.65 | - | 1.01/0.92 | 0.90/0.82

Interpretation:
- POSTREFIT-EXACT-H and PROFILE-DQR produce nearly identical estimates
  (A_k same order), so residual global Dantzig regularisation is NOT the
  blocker; the exact-row sandwich is merely noisier (se/sd > 1.1).
- POP-H is in band at n=400 (0.94-0.98) and near-band at n=200 (0.88-0.96):
  the h_inf / population-direction layer is healthy.
- The remaining practical gap is the FITTED-direction sandwich variance:
  se/sd 0.27-0.68 (n=200) -> 0.62-0.80 (n=400), improving with n; coverage
  0.72-0.86 at n=400 vs band [0.88,0.98].

## 3. Decisions per actions

- Item 16: beta_refit error is now small (0.25-1.6) and A_k is small
  (0.2-1.2) => the mu rate is NOT reopened.
- Item 17: P06-v5 moves materially toward POP-H/TRUE-SUPPORT (A_k ~0.2,
  coverage 0.72-0.86) but remains outside the gate; P05-to-P06 scaling is
  coherent (l1_refit 1.57 -> 0.25, se/sd 0.27-0.68 -> 0.62-0.80).
  => P07-v5 (n=800, B=50) is conditionally authorised as the scaling
  locator; launched under this report. P07-B=200 still requires review.
- Item 18: the full P01-P06 rerun remains unauthorised pending review of
  this report.

## 4. Round-6 candidate question (not a ruling)

The remaining under-coverage is a finite-sample variance-layer gap of the
practical fitted-direction sandwich (se/sd 0.62-0.80 at n=400) that POP-H
(0.90-1.09) does not show. Whether the manuscript variance statement needs a
theory-sanctioned amendment (e.g., a stated higher-order term in the fitted
meat) is deferred until P07-v5 scaling is measured.

## 6. P07-v5 scaling locator (n=800, B=50) — addendum

Run: results/raw/simulation/pilot_v5_P07_B50 (4 shards, 0 runner errors,
150 method rows; methods PROFILE-DQR | POSTREFIT-EXACT-H |
PROFILE-DQR-TRUE-SUPPORT; no POP-H at n=800 because no population-direction
asset exists for the P07 cell).

First stage: l1_lasso med 1.569 -> l1_refit med 0.163 (q25 0.135, q75 0.209),
l2_refit 0.075, sel_size med 5 (all actives), refit_set 7, omitted 0,
post_refit ok 1.000, status ok 1.000.

PROFILE-DQR: bias <= 0.003; coverage 0.840/0.840/0.840/0.940 (mcse 0.03-0.05);
se/sd 0.77/0.72/0.77/0.87. TRUE-SUPPORT: coverage 0.88-0.96, se/sd 0.89-0.99.

Three-cell nested scaling of the practical variance gap (active coords):

| cell | se/sd x00001 | se/sd x00002 | cov x00001 | cov x00002 |
|---|---|---|---|---|
| n=200 | 0.68 | 0.27 | 0.820 | 0.360 |
| n=400 | 0.75 | 0.65 | 0.860 | 0.760 |
| n=800 | 0.77 | 0.72 | 0.840 | 0.840 |

Interpretation: the fitted-direction sandwich under-estimation closes
monotonically with n (se/sd 0.27-0.68 -> 0.72-0.87) but remains ~13-28% below
1 at n=800, keeping practical coverage at 0.84-0.94 (1 mcse below the lower
band edge for 3/4 coordinates) while TRUE-SUPPORT is in band. Bias and
precision are no longer limiting at any cell.

## 7. Status and next decision

Per action 19: P05 fails while P06/P07 move toward/into calibration but do not
fully enter the band at B=50. The full P01-P06 rerun and P07-B=200 remain
unauthorised pending review of this report. Candidate round-6 questions:
(i) whether the persistent ~15-25% fitted-meat shortfall warrants a
theory-sanctioned higher-order variance statement (POP-H and TRUE-SUPPORT are
in band, localising the term to the fitted-direction meat);
(ii) the manuscript inference-regime wording (n>=400/800) and the GSE65391
confirmatory high-dimensional Wald scope, per round-4 Q5.

## 5. Compliance

No mu/lambda/h change; round-5 spec implemented verbatim; all 300 method
rows per cell accounted for; 0 runner errors; truth used only in the
registered simulation diagnostics (l1/l2 errors, omitted-active count,
A_k_lasso_start).

Evidence: results/raw/simulation/pilot_v5_P05_B50/ and pilot_v5_P06_B50/
(replication_metrics, coordinate_metrics, theory_diagnostics,
screening_records + merge_manifest.json).
