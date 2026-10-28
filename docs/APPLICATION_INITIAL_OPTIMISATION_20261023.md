# Initial L1 numerical optimisation

The prediction and independent gene-selection procedures use the complete convex
profile criterion at `h_est`, with the frozen two-pass cluster-score loadings and
the unchanged coordinate penalty. The numerical option
`JPDMEQR_INITIAL_ACCELERATION=clinical_kkt` checks a restricted solution using the
entire fitted design. Its alternative value `none` retains the proximal path.

Let `U` contain the coordinates whose penalty factor equals zero and `P` contain
the remaining coordinates. Each loading pass solves the complete unpenalised
profile criterion in `X_U` using `fit_profile_unpenalised_newton_v2`, sets the
coordinates in `P` to zero, and evaluates the original full-design profile
gradient. The clinical Newton KKT tolerance is the smaller of `1e-9` and the
caller tolerance. Its coefficient-change tolerance is the smaller of `1e-9`
and the caller tolerance. The existing strict nuisance tolerance and profile
identity check remain part of acceptance.

The full-design certificate retains its score, cluster-score matrix, complete
Hessian and identity diagnostics. It does not request the four fold-specific
Hessian sums, which are used by precision-row cross-validation. The prediction
initial estimator has no precision-row stage; these fold fields are `NULL`.

For a coordinate penalty `a_j`, the full-design certificate evaluates

\[
r_j=|g_j|\quad(j\in U),\qquad
r_j=(|g_j|-a_j)_+\quad(j\in P).
\]

Both the maximum original-unit residual and the existing normalised KKT
residual must be at most `min(caller tolerance, 1e-7)`. Convexity makes the
zero-residual condition sufficient for global optimality of the full L1
criterion. The recorded positive tolerance gives the declared numerical KKT
accuracy. When a gene coordinate exceeds its penalty bound, the original
full-design proximal solver continues from the computed clinical coefficients.
All original solver failures terminate their request and retain the condition.

Pass-zero loadings are always computed with every coefficient equal to zero.
The second loading vector is computed at the preliminary fit, and the final
coordinate loadings are their coordinatewise maximum. Each optimisation uses
its fixed loading vector. Neither the bandwidth nor either penalty component
depends on a validation outcome or on the numerical certificate.

The clinical request may be reused within one training fold. Cache equality
requires identical observed outcomes, clinical design, nuisance design, ordered
subject identifiers, quantile, bandwidth, `Lambda`, numerical controls and solver
configuration. `12_validate_GSE65391.R` creates a separate environment for each
inner training fold. The request does not contain a validation response.

Every fit retains its acceleration setting, full coordinate residuals, clinical
Newton history, proximal history when executed, and cache-use indicators.
`scripts/research/check_application_initial_clinical_20261023.R` saves actual
reference and accelerated fits, their coefficient and prediction differences,
selection indicators, loading vectors and elapsed times. The complete timing
grid contains all 108 frozen candidates on the recorded observed training fold.

The matched low-dimensional workflow uses zero L1 penalties throughout and
does not call these initial-estimator or prediction-candidate functions. Its
source comparison excludes only the named top-level declarations in application
parts 10 and 13. The remaining top-level expressions are compared in full using
the R parser. Each matched job retains the names and reasons for these exclusions.
