# PROFILE-DQR-SPLIT

This procedure screens on one deterministic cluster subset and fits the
profile estimator on its disjoint complement. Quantile-score screening uses
only screening responses. Prespecified inference coordinates are retained
according to the declared forced-inclusion rule. The selected design is
then fitted with the same complete profile, paired score, paired precision
validation and cluster sandwich as PROFILE-DQR.

For an already prepared split, the adapter requires disjoint role sets and
checks that every actual training cluster belongs to the inference role.
It also checks that the training observations contain no screening cluster.
Returned coefficients and predictions use the original named columns; the
screening record retains both role sets and the selected global coordinates.

Allowed metrics are estimation, selection, prediction and coordinate coverage
for the independent inference sample. The paired benchmark test reconstructs
that exact sample and design and checks agreement with a direct fit. It
verifies deterministic roles and rejects full-data input with a split-role
record. Dependencies and the acceptance entry point follow the PROFILE-DQR
source note.
