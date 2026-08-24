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

## 5. Compliance

No mu/lambda/h change; round-5 spec implemented verbatim; all 300 method
rows per cell accounted for; 0 runner errors; truth used only in the
registered simulation diagnostics (l1/l2 errors, omitted-active count,
A_k_lasso_start).

Evidence: results/raw/simulation/pilot_v5_P05_B50/ and pilot_v5_P06_B50/
(replication_metrics, coordinate_metrics, theory_diagnostics,
screening_records + merge_manifest.json).
