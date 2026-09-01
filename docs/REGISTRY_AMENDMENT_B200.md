# Registry amendment: main registry at B=200 (versioned)

**Status:** project-level decision by the study lead (2026-09-01), recorded
prospectively; new versioned registry `config/simulation_main_b200.csv`;
frozen file `config/simulation_main.csv` unchanged. Escalation pathway per
AGENTS.md §8.3 ("a necessary change creates a new versioned registry and run
ID").

## Decision

Run the final main registry at B=200 for all inferential rows (modules A-H
excluding computation module G). Module G rows keep their frozen B=50
(computation-scaling cells). All 78 rows, experiment IDs, methods, seeds,
tuning grids, targets and DGP branches are otherwise unchanged.

## Reason

Compute/time budget: the frozen B=500/1000 scale evaluates to ~95,000
core-hours (6-9 weeks at 64-core full load); B=200 evaluates to ~29,500
core-hours (~10-12 days), freeing time and budget for the application layer
and manuscript completion. No statistical object, method, threshold, or
design branch changes.

## Disclosure and precision

- Every reported Monte Carlo mean/proportion continues to carry MCSE
  (binomial intervals for coverage/type-I), per the results contract.
- At B=200, coverage MCSE ~ sqrt(0.25/200) = 0.035 (worst case); the pilot
  (also B=200) already established the same precision for the calibration
  tables.
- The manuscript will state "B=200 Monte Carlo replications throughout
  (computation cells B=50)".
- The pilot record (P01-P06, B=200) already provides identical-scale
  mechanism evidence; nothing in the freeze evidence changes.

## Consistency with governance

- No threshold change; no method/target/DGP change; no comparator change.
- Run ID namespace: final runs use run_id with a `_b200` suffix to keep
  provenance distinct from any B=500/1000 run.
- Theory-side visibility: this amendment is committed to the branch for
  GPT-5.6 Pro review; no theory statement depends on the Monte Carlo scale
  (asymptotic statements are n-asymptotic, not B-asymptotic).