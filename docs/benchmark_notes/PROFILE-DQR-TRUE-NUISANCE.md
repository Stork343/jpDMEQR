# PROFILE-DQR-TRUE-NUISANCE

This diagnostic addresses the structural coefficient when the generated
nuisance design is correctly specified. It subtracts the actual generated
`Z_i b_i` from the response, matches random-effect rows by cluster identifier,
and fits the resulting fixed-effect problem. A zero nuisance design makes
its nuisance optimum exactly zero. The fixed-effect score and curvature are
therefore the ordinary smoothed quantities averaged within each cluster.
The paired combination, Dantzig rule and cluster sandwich then apply to
that known-nuisance problem.

The paired adapter rejects a mismatch between generated nuisance dimension
and fitted design, and rejects an additional nonlinear generated component.
Simulation scenarios using this diagnostic must declare a correctly specified
structural target. The practical adapter receives no generated random effects.

Allowed metrics are estimation, prediction, structural-coordinate coverage
and the nuisance-estimation gap. The paired benchmark test independently
subtracts the generated nuisance and compares the direct fitted problem,
including score and variance calculations. It also tests dimension rejection.
Implementation, dependencies and acceptance entry point follow the
PROFILE-DQR source note.
