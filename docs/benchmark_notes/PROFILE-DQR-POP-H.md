# PROFILE-DQR-POP-H

This diagnostic retains the practical selection, reduced refit and empirical
paired score. It replaces the estimated precision direction by a direction
from an independent population sample. The population Jacobian is the paired
complete-profile Jacobian at the analysis bandwidth, evaluated at the
declared scenario target. Its covariance combines the two bandwidth scores
on each population cluster before centring.

The adapter requires the scenario configuration through `context$config`,
the actual implementation commit through `control$implementation_commit`,
and an explicit final-asset policy. The validator checks the scenario,
full-design coordinate coverage, bandwidths, working Lambda, construction
version and implementation identity. Final assets require at least 100000
clusters and four independent repetitions together with the declared
convergence and population diagnostics.

Allowed metrics are coordinate coverage and the population-precision gap.
The paired benchmark test constructs genuine development target and direction
assets from independent generated samples. It checks the population inverse,
direct score and sandwich calculation and rejection of mismatched scenario,
commit and asset scale. This test certifies the adapter calculation; formal
population objects pass their own final-scale gate. Dependencies and the
acceptance entry point follow the PROFILE-DQR source note.
