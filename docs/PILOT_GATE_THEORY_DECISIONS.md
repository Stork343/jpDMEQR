# Pilot-gate theory decisions

The authoritative detailed decision is stored in `results/preflight/PILOT_GATE_THEORY_DECISION.md` at the same commit. This short document promotes its status into the documentation hierarchy so subsequent agents do not treat the failed-pilot questions as unresolved.

Until `docs/METHOD_SPECIFICATION.md`, code, configs and tests are synchronised, the following decisions supersede the old provisions on sandwich meat, primary inference bandwidth and Dantzig tuning:

1. Primary Wald meat for the **unsmoothed** regularised profile target uses the ordinary quantile score
   `psi_tau(u)=tau-I(u<0)` at fitted profile residuals. The smoothed-score meat is retained only as a diagnostic/sensitivity. The bread remains the corrected smoothed profile Hessian.
2. Replace the primary inferential bandwidth `h=c_h n^{-1/3}` by `h=c_h n^{-3/10}`, `c_h in {0.75,1,1.25}`. This gives both `sqrt(n) h^2 -> 0` and `n h^3 -> infinity`.
3. Expand the Dantzig multiplier candidates to `c_mu in {0.10,0.25,0.50,1,2,4}`. Feasibility alone is insufficient. Select the constant without truth/coverage using the frozen cluster-level inverse-Hessian cross-validation rule in the detailed decision and add inverse-defect/POP-H row diagnostics.
4. No scalar cluster-size correction and no empirical SE multiplier are permitted.
5. Existing pilot thresholds remain unchanged until the corrected pipeline is rerun.
6. A pure `mu` change does not mathematically change target/direction assets, but this decision also changes `h`; therefore POP-H directions must be rebuilt. Profile-target assets may be revalidated if their mathematical dependency hash is unchanged.
7. Fix POP-H construction so its analysis bandwidth is determined by the **scenario analysis cluster count**, not by the population Monte Carlo sample size.
8. Retain Module-A `n=100` cells as stress/scaling cells; do not claim near-nominal inference unless the corrected results support it.

No final-scale simulation or confirmatory empirical inference is authorised until the detailed implementation sequence in `results/preflight/PILOT_GATE_THEORY_DECISION.md` is completed and a new pilot gate passes.