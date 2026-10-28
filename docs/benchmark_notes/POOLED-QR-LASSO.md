# POOLED-QR-LASSO

The adapter solves the pooled check-loss sum plus the declared L1 penalty
on all supplied fixed-effect coefficients. It uses quantreg's Frisch Newton
linear-programme solver and the exact opposite-diagonal representation of
the absolute-value terms. The numerical solution and coefficient names are
retained. The package documentation is available at
<https://CRAN.R-project.org/package=quantreg>.

The reference target is a pooled conditional quantile coefficient. The
adapter ignores the cluster nuisance effects and returns no interval or
cluster sandwich. The runner chooses its penalty using training clusters;
new-cluster predictions use the named fixed-effect coefficient vector.
Allowed metrics are estimation, selection and prediction.

The acceptance tests compare the adapter with an independently expressed
CVXR linear programme on identical observations and coordinate penalties.
They check deterministic predictions and cluster-label invariance and verify
that inferential fields are absent. The quantreg package has GPL licensing. Current acceptance uses
the actual installed version recorded in the session information.
