# POSTREFIT-EXACT-H

This simulation mechanism diagnostic retains the practical convex selection
and reduced refit. On the selected refit coordinates it computes the paired
Jacobian inverse exactly, embeds its requested columns in the full design,
and uses those directions with the paired score and its cluster sandwich.
The selected paired Jacobian must be positive definite. The direction has
zero coordinates outside the reduced design.

The construction measures the effect of global precision estimation relative
to inversion on the refit design. It is an oracle mechanism comparison.
Its declared metrics are estimation and the precision gap. The paired
benchmark test verifies the inverse equations on the actual fitted Jacobian,
the zero entries outside the reduced design, deterministic output and the
same-score correction and variance identities. The formulas, dependencies
and acceptance entry point follow the PROFILE-DQR source note.
