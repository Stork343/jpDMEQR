# Pilot-gate theory decisions

The current authoritative decision is:

`results/preflight/PILOT_V5_THEORY_DECISION_ROUND6.md`

with the executable specification in:

`docs/METHOD_SPECIFICATION_ROUND6_AMENDMENT.md`.

Round 6 supersedes only the variance-gate interpretation, pilot authorisation path and empirical inferential scope. The point estimator and all non-conflicting round-1 through round-5 method rules remain in force.

Authoritative rules:

1. The round-five point estimator is unchanged: score-loaded L1 selection at `h_est`, zero-L1 post-profile refit at `h_inf`, full-design reprofile, frozen Dantzig selector, one-step update, and unsmoothed cluster-score primary meat.
2. The P07-v5 `B=50` result (`SE/SD≈0.72--0.87`, coverage `0.84--0.94`) is **not** declared calibrated. It is a material finite-sample discrepancy, but it does not identify a universal missing variance multiplier.
3. No scalar SE inflation, `mu`-slack variance addition, generic KC/CR2/MD/CR3 production correction, or selected-model df factor is authorised.
4. Before any full pilot, run the required P06-v5 `B=50` `T=L+Q+R` decomposition using the existing POP-H asset:
   - `L` = population-direction first-order influence term;
   - `Q` = finite-sample direction-estimation contribution `(omega_hat-omega_pop)' G_star`;
   - `R` = all remaining one-step remainder.
   Verify the complete variance/covariance identity.
5. Compute four variance-ladder levels: fitted direction/fitted score, fitted direction/target score, population direction/fitted score, population direction/target score. Compare them to empirical `Var(L)`, `Var(L+Q)` and `Var(T)`.
6. For mechanism ratios use 5,000 paired bootstrap resamples of the replication index, seed `20260825`, to quantify Monte Carlo uncertainty. Bootstrap is diagnostic only and never changes a model confidence interval.
7. Mechanism classification:
   - V1: fitted meat matches `Var(L+Q)` but `Var(T)` is larger because of `R`/covariance -> keep current asymptotic sandwich and report finite-sample remainder;
   - V2: fitted meat underestimates `Var(L+Q)` while population-direction meat is calibrated -> direction-estimation variance is missing; stop for a direction-aware variance amendment;
   - V3: population-direction fitted-score meat is low while target-score meat is calibrated -> fitted-score/profile plug-in squeeze; stop for a profile-specific correction derivation;
   - V4: mixed/inconclusive -> add 50 P06-v5 diagnostic replications only and repeat the ladder.
8. Full P01--P06 and P07-B=200 are blocked until the variance mechanism is classified. Under V1 only, P07-B=200 is authorised first; if no new anomaly appears, the full versioned P01--P06 B=200 pilot is then authorised automatically under the same variance rule.
9. Under V1, the SE/SD and coverage bands remain visible finite-sample calibration diagnostics but cease to be binary implementation-correctness blockers for intentionally small-n stress cells. Under V2--V4 the variance mechanism remains unresolved and the full pilot stays blocked.
10. GSE65391 at approximately 129 independent subjects does **not** support confirmatory high-dimensional gene-level Wald inference under current evidence. High-dimensional analysis is restricted to prediction, stability and exploratory rankings/intervals.
11. A prespecified low-dimensional application layer is allowed in principle: direct unpenalised profile fit of frozen clinical covariates/immune modules with `d_app<=floor(n/log n)`, reporting ordinary cluster-sandwich and delete-one-subject jackknife intervals side by side. Claim strength is conditional on application-matched simulations.
12. Authorise PAPP-HD (`n=129,p=500,s=5,tau=.5,B=200`, empirical GSE visit-count design) and PAPP-LD25/50/75 (`n=129,d=20,s=5,B=200` each, direct low-dimensional profile fit) as application-scope diagnostics. No subject identifiers or outcomes enter the DGP; only the frozen aggregate visit-count design is used.
13. P05/P06-v4 B=50, P05/P06-v5 B=50, P07-v5 B=50 and the round-six ladder must appear in the final pilot manifest under `mechanism_evidence` with `formal_gate=false`.
14. Final execution order is variance classification -> any required method amendment -> final code/config freeze -> tests -> strict geometry -> micro-preflight -> authorised formal pilot/scaling runs -> pilot manifest -> freeze-preflight -> B=500/1000 main registry.
15. Final-scale simulation and confirmatory empirical inference remain blocked until those conditions are met.
