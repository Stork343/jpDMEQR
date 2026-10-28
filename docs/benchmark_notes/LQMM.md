# LQMM

The comparator directly uses the maintained lqmm package for linear quantile
mixed models. Its reference is Geraci and Bottai (2014), Linear quantile mixed
models, Statistics and Computing 24, 461–479,
DOI `10.1007/s11222-013-9381-9`. Package documentation is available at
<https://CRAN.R-project.org/package=lqmm>.

The adapter supplies a fixed-effect formula with the named design columns
and a random intercept or intercept with time slope, matching its declared
working nuisance design. It uses Gaussian quadrature and the recorded
control values. Its design must be a declared common low-dimensional set
or an explicitly oracle support. Prediction for a new cluster uses the
fixed-effect coefficients.

The current adapter returns point estimates and a serialisable lqmm fit.
Its allowed outputs are low-dimensional estimation and prediction; it does
not construct a coordinate interval. Practical independently screened and
oracle-support uses have separate labels. The lqmm package uses GPL licensing.
Fidelity checks compare the adapter with a direct package call on the same
formula, observations and controls, and examine a negligible random-effect
variance limit. Current package version and numerical conditions are recorded
with acceptance results.
