# Pilot-gate theory decisions

The current authoritative decision is:

`results/preflight/PILOT_V2_THEORY_DECISION_ROUND3.md`

with the executable specification in:

`docs/METHOD_SPECIFICATION_ROUND3_AMENDMENT.md`.

Round 3 supersedes the first-stage penalty rule and the conflicting action sequence. Non-conflicting round-1 and round-2 decisions remain in force.

The following rules are authoritative:

1. The inferential target remains the **unsmoothed regularised profile parameter**. The exact profile score, corrected Schur-complement Hessian, one-step formula and POP-H definition are unchanged.
2. The raw unit-loaded rule `lambda_beta=sqrt(log p/n)` is revoked. The first stage uses coordinate penalties calibrated to the empirical standard deviation of the exact smoothed **independent-cluster profile score**.
3. Freeze
   \[
   \alpha_{\lambda,n}=0.10/\log\{\max(n,3)\},
   \quad
   q_{\lambda,n}=\Phi^{-1}\{1-\alpha_{\lambda,n}/(2p_P)\},
   \quad
   \lambda_{0,n}=1.10q_{\lambda,n}/\sqrt n.
   \]
   At `b=0`, compute centred cluster-score loadings `ell_j^(0)`, fit a preliminary estimator, update once at `b^(1)`, set `ell_j^final=max(ell_j^(0),ell_j^(1))`, and refit with coordinate penalty `lambda_0,n ell_j^final`. No lambda CV and no truth/coverage tuning are allowed.
4. Implement the loaded penalty through scalar `lambda_beta=lambda_0,n` and vector `penalty_factor_j=ell_j^final*base_penalty_factor_j`. Non-finite or numerically degenerate retained loadings cause an explicit failure.
5. A final first stage is accepted only when the weighted normalised KKT residual is at most `1e-3`, maximum nuisance gradient is at most `1e-7`, and the last coefficient change is at most `1e-7*max(1,||beta_hat||_inf)`. Use at least 2000 outer iterations and 50 backtracking steps; non-convergence remains a failed replication.
6. The sparse first-stage rate remains `||beta_hat-beta_h_star||_1=O_p{s sqrt(log p/n)}` under score domination, uniformly comparable loadings, profile RSC and numerical KKT error `o_p(lambda_0,n)`. For the unsmoothed target add the smoothing-target discrepancy. The rate is not a unit-constant finite-sample cutoff, but an order-one empty fit is outside the intended regime.
7. `lambda_beta_multipliers` now denotes explicit sensitivity multipliers around the calibrated coordinate penalty. The primary method uses multiplier `1`; alternative values occur only in separately registered sensitivity cells.
8. Round 3 makes **no change** to the Dantzig anchor/grid/selector. Keep `h=c_h n^{-3/10}`, `c_mu in {0.02,0.05,0.10,0.25,0.50,1,2,4}`, two-fold cluster CV for `n<200`, four folds for `n>=200`, held-out inverse-defect loss, and smaller-`mu` tie-breaking.
9. The precision stage begins only after the final lambda-calibrated fit passes. Recompute nuisance profiles and all Hessian contributions at that final fit before selecting `mu`. Reopen the `mu` rate only if actual `A_k`, POP-H row errors, total Bahadur remainder and P05--P06 scaling remain poor jointly after lambda repair; `D_k` alone is not a hard gate.
10. The primary variance remains corrected smoothed bread plus unsmoothed fitted cluster-score meat. The P05 TRUE-SUPPORT result (`SE/SD=0.80--0.95`, coverage `0.87--0.92`, negligible bias) is acceptable as mild finite-sample diagnostic behaviour and does not authorise a scalar correction or reopen the variance formula. Jackknife/KC/CR2/MD/CR3 remain diagnostics.
11. The practical proposed-method P05 gate retains `SE/SD in [0.80,1.20]` and coverage in `[0.88,0.98]`; P06 confirms scaling. Oracle methods are mechanism diagnostics and do not independently fail the freeze because one oracle coordinate is marginally below the coverage band.
12. The lambda change invalidates all practical P01--P06 results for freeze purposes. Create a versioned pilot registry and new run IDs, run a P05 `B=50` lambda-mechanism diagnostic first, then rerun P01--P06 from scratch if the first stage is healthy.
13. Profile-target and POP-H assets may be reused under matching dependency hashes because analysis lambda is not a mathematical dependency. Practical fits, Hessians, precision rows, intervals and summaries must be regenerated.
14. No final-scale `B=500/1000` simulation or confirmatory empirical inference is authorised until the round-three action sequence is completed and a fresh matching freeze manifest passes.
