# Pilot-gate theory decisions

The current authoritative decision is:

`results/preflight/PILOT_V4_THEORY_DECISION_ROUND5.md`

with the executable specification in:

`docs/METHOD_SPECIFICATION_ROUND5_AMENDMENT.md`.

Round 5 supersedes the direct use of the penalised L1 coefficient vector as the one-step starting estimator. Non-conflicting round-1 through round-4 decisions remain in force.

Authoritative rules:

1. The round-four paired probe confirms that the dual-bandwidth first stage is working: at `n=400` support recovery is exact and RSC diagnostics are healthy, but the penalised `l1_error` remains about `2.27`. Exact support therefore does not imply a sufficiently small starting-estimator error.
2. The round-four precision-review condition is met, but a faster asymptotic `mu` rate is **not authorised**. The leading `sqrt(log p/(n h_inf))` term is the kernel-Hessian max-norm fluctuation scale; existing debiased-SQR theory uses a sufficient tolerance containing this term plus additional plug-in/nonlinear terms rather than a faster rate.
3. Keep the round-two Dantzig anchor, candidate grid and held-out inverse-defect CV selector unchanged during the round-five mechanism probe.
4. Add a **post-L1 unpenalised profile refit**. The round-four penalised estimator remains the selection estimator. Define the refit set as the union of its nonzero penalised coordinates, all always-included coordinates and all registry-prespecified inference target coordinates.
5. Require refit dimension `d_R <= floor(n/log(max(n,3)))`. Failure is explicit; there is no fallback to the penalised start.
6. Refit with zero L1 penalty at `h_inf`, using the same nuisance ridge. Require reduced-model gradient `<=1e-7`, nuisance gradient `<=1e-7`, coefficient change `<=1e-8*max(1,||beta_refit||_inf)`, positive reduced Hessian minimum eigenvalue `>1e-8`, and condition number `<1e10`.
7. Embed `beta_refit` into the full `p` vector, then recompute every full-design nuisance profile, exact score, effective Hessian and fold Hessian at `h_inf`. Only these post-refit full-p objects may enter Dantzig and one-step inference.
8. The practical one-step is now `beta_refit_k - omega_hat_k' g_hinf(beta_refit)`. The original penalised vector is retained separately as `beta_l1` for selection and diagnostics.
9. Under sure inclusion and `d_R=O(s)`, the post-refit has schematic low-dimensional error `||beta_refit-beta_star||_1 = O_p{d_R/sqrt(n)+d_R h_inf^2}`. The normalized inverse-defect bound then loses the extra `sqrt(log p)` factor created by the penalised-start L1 error. If active variables are omitted, an omitted-signal approximation term remains.
10. Add the simulation-only `POSTREFIT-EXACT-H` diagnostic: invert the data-selected reduced Hessian on the refit set, extend the requested row by zero, and form the one-step coordinate from the full-p score. It may not tune `mu` or appear as a primary method.
11. The observed `x00002` bias at exact support is not a new first-order influence term. It is consistent with the non-negligible remainder `(e-H omega)' Delta` because the penalised `Delta` remains large. The post-refit probe is the required discriminator.
12. Do not run the full P01--P06 under the current penalised-start rule and do not run P07-B=200 now. After implementing/tests, run paired P05-v5/P06-v5 at `B=50` with nested clusters/common random numbers where feasible.
13. If P05/P06-v5 substantially reduce refit error and `A_k`, and P06 moves toward POP-H/TRUE-SUPPORT but remains outside the gate, P07-v5 is conditionally authorised first at `B=50`, not `B=200`.
14. The practical calibration bands remain `SE/SD in [0.80,1.20]` and coverage in `[0.88,0.98]`. They are freeze/debug criteria, not theorem promises. No formal gate is issued from `B=50`.
15. Profile-target and POP-H assets remain reusable under matching dependency hashes because target, DGP and `h_inf` are unchanged. All practical refits, full inferential reprofiles, precision rows, intervals and summaries must be regenerated.
16. Final `B=500/1000` simulation and confirmatory GSE65391 high-dimensional inference remain blocked.
