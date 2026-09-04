# P01 mechanism check (B=50) — round-2 decomposition update

Run: results/raw/simulation/pilot_v2r2_P01_B50b (P01, B=50, dev settings)
Pipeline: round-2 (inverse-defect mu selector, 9-candidate grid, A_k, D_k demoted,
variance ladder, CR0/CR3/KC FS diagnostics), commit 399f17b.

## 1. The new mu selector did NOT close the precision defect at n=80
- Selected mu: c=0.5 (mu=0.284) in 109/200 rows, c=1.0 (mu=0.569) in 91/200.
  The defect-loss CV still lands on c>=0.5; smaller candidates (0.02/0.05/0.10)
  are infeasible on the noisy 2-fold training Hessians at n=80, so they are
  excluded and the argmin falls back to 0.5-1.0.
- A_k (actual normalized inverse-defect inner product): median 3.09, q90 6.99.
  D_k (Holder bound): median 11.39. E2: median 0.59. cosine 0.89.
  => The precision rows genuinely corrupt the first-order expansion (A_k ~ 3-7,
  not just a loose Holder bound).
- PROFILE-DQR P01 (B=50): bias -0.40/+0.64, coverage 0.00-0.32, se_sd_ratio
  0.48-0.88. TRUE-SUPPORT: bias ~0.01, coverage 0.74-0.86, se_sd_ratio 0.60-0.87.

## 2. The variance ladder localises the oracle SE shortfall
For PROFILE-DQR coordinates (n=50 reps):
  n*se_emp^2 ~= Var(T) (=0.36-0.40) by construction.
  Var(L) (oracle influence at target) ~ 0.56-0.74  (> Var(T)).
  Var(R) (remainder) ~ 0.53-0.67   (first-order magnitude, NOT negligible).
  2Cov(L,R) ~ -0.74 to -1.02 (strongly negative; L and R cancel in T=L+R).
  ladder se1 (population) ~= se2 (target-score) ~= se3 (fitted-score pop dir)
  ~= 0.104-0.105 — the meat level is stable across population/target/fitted.
  se4 (fitted + sample exact inverse) NaN (exact solve not computed in this band).

Interpretation: the sandwich underestimates Var(T) because (a) the selected
(imprecise) omega shrinks the projected scores relative to the population
direction, and (b) Var(R) is of the same order as Var(L), so the one-step
estimator is NOT dominated by the oracle influence at n=80; the CR0/CR3/KC
leverage corrections only inflate the sandwich by 3-14%, insufficient for the
gap. This confirms: the deficit is not a classic small-sample sandwich
downward bias; it is the combination of imprecise precision rows at n=80 and a
non-negligible remainder term.

## 3. FS diagnostics (TRUE-SUPPORT)
  CR0->CR3 inflation only ~14% (0.045->0.051), ~3% for KC. Leverage-type
  corrections do not bridge the oracle gap.

## 4. Mechanism status
The diagnostic pipeline is coherent (all components run, T=L+R identity holds,
ladder levels well-behaved), but the mechanism is NOT closed at n=80: precision
rows are non-first-order (A_k~3-7) and the remainder is large. The decisive
calibration cell is P05 (n=200); the full P01-P06 run (frozen replication
counts) is required before any gate decision. No threshold or correction has
been altered.

Evidence: results/raw/simulation/pilot_v2r2_P01_B50b/{coordinate_metrics,
replication_metrics,theory_diagnostics}.csv + summaries.
