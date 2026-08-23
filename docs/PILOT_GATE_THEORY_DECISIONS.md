# Pilot-gate theory decisions

The current authoritative decision is:

`results/preflight/PILOT_V3_THEORY_DECISION_ROUND4.md`

with the executable specification in:

`docs/METHOD_SPECIFICATION_ROUND4_AMENDMENT.md`.

Round 4 supersedes the use of the inference bandwidth inside the penalised first stage. Non-conflicting round-1, round-2 and round-3 decisions remain in force.

Authoritative rules:

1. The P05/P06 baseline DGP is unchanged. No reduction of AR(1) correlation, signal change, cluster-size change, or nuisance removal is authorised. The current evidence does not reject population identification.
2. The single-bandwidth first-stage RSC guarantee is not certifiable at `n=200` or `n=400`. The standard high-dimensional SQR sufficient lower order contains `sqrt{s log(p)/n}`; for P05 this is about `0.394` versus the inference bandwidth `0.204`, and for P06 about `0.279` versus `0.166`.
3. Introduce a **dual-bandwidth** procedure. The first-stage estimation bandwidth is
   \[
   h_{est}=\{\log(p_P\vee2)/n\}^{1/4},
   \]
   where `p_P` is the number of penalised fitted coordinates. The primary constant is `c_est=1` and is not data-selected.
4. Keep the inference bandwidth unchanged:
   \[
   h_{inf}=c_hn^{-3/10},\qquad c_h=1
   \]
   for the primary pilot.
5. Compute the round-3 cluster score loadings and both penalised first-stage passes using `h_est`. Keep the round-3 lambda critical value, two-pass loadings, weighted KKT contract, and primary lambda sensitivity multiplier `a=1` unchanged.
6. Do **not** lower `a` or add lambda CV/BIC in round 4. The P05 evidence was obtained in an uncertified first-stage curvature regime; changing lambda before repairing that regime would confound mechanisms and weaken the retained score-domination argument.
7. After the final first stage passes, discard its nuisance/Hessian objects for inference. At the accepted `beta_hat`, reprofile nuisance effects and recompute the exact score/Hessian/fold Hessians using `h_inf`; only then run the frozen Dantzig selector and one-step correction.
8. The first-stage bound is written
   \[
   ||beta_hat-beta_star||_1
   =O_p\{s sqrt(log p/n)+s h_{est}^2\}.
   \]
   Since `h_est^2=sqrt(log p/n)`, the retained order is still `O_p{s sqrt(log p/n)}` up to constants.
9. The Dantzig anchor/grid/selector is unchanged from round 2. Reopen `mu` only if, after the dual-bandwidth first stage is healthy, actual `A_k`, POP-H row errors, total Bahadur remainder and P05-to-P06 scaling remain jointly poor. `D_k` alone is not a hard gate.
10. The primary variance remains corrected smoothed inferential bread plus unsmoothed fitted cluster-score meat. No scalar finite-sample variance correction is authorised.
11. Add simulation-only first-stage RSC diagnostics: population/sample active-set minimum Hessian eigenvalues, condition numbers, a deterministic cone-curvature proxy, and `h_est/sqrt{s log(p)/n}`. Truth is diagnostic only and cannot tune the method.
12. Create a versioned `pilot_v4` registry. After implementation/tests, run **both** P05-v4 and P06-v4 at `B=50` as a paired mechanism/scaling probe. The old P05-v3 result is historical only because the estimator has changed.
13. The B=50 probe is not a formal coverage gate. Do not launch the full P01--P06 rerun until its first-stage/RSC/Bahadur scaling is reviewed.
14. The full practical P05 gate remains `SE/SD in [0.80,1.20]` and coverage in `[0.88,0.98]`; these are freeze/debug criteria, not theorem promises.
15. If P05/P06 show coherent monotone scaling but the later full P05 gate still fails, a P07 scaling cell is conditionally authorised: `n=800,p=500,s=5,tau=0.50,q=1`, same baseline DGP, `B=200`. P07 does not automatically replace P05 as the gate; any change in manuscript applicability requires a later theory/project decision.
16. Profile-target assets remain reusable under matching dependency hashes. POP-H remains defined by `h_inf` and is reusable if its inferential dependency hash is unchanged. Every practical fit/Hessian/direction/interval must be regenerated.
17. No final-scale `B=500/1000` simulation or confirmatory GSE65391 high-dimensional inference is authorised until the round-4 mechanism probe and subsequent freeze sequence pass.
