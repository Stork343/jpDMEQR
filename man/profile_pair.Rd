\name{profile_pair}
\alias{profile_pair}
\alias{fit_profile_pair_v2}
\alias{infer_profile_pair_v2}
\alias{profile_bandwidth_pair_v2}
\alias{build_profile_pair_target_asset_v2}
\alias{build_profile_pair_direction_asset_v2}
\alias{validate_profile_pair_target_asset_v2}
\alias{validate_profile_pair_direction_asset_v2}
\alias{profile_pair_target_scores_v2}
\title{Paired bandwidth profile inference and population assets}
\description{
Fits regularised mixed-effects quantile profiles and constructs coordinate
intervals from paired bandwidth scores. Population builders provide numerical
profile targets and matched precision directions for registered simulation
designs. The method identity is \code{profile_pair_20261022}.
}
\usage{
fit_profile_pair_v2(y, X, Z, cluster_id, tau, h_inf,
  h_est = NULL, lambda_gamma = 1, mean_score_multiplier = 1,
  lambda_multiplier = 1, penalty_factor = NULL,
  target_coordinates = NULL, control = list())

infer_profile_pair_v2(fit, coordinates, mu_grid, ci_level = 0.95,
  solver_preference = c("CLARABEL"), solver_opts = list(),
  df_correction = FALSE)

profile_bandwidth_pair_v2(y, X, Z, cluster_id, beta, tau, h,
  lambda_gamma = 1, nuisance_control = list(),
  identity_check = FALSE, assembly = NULL, need_fold_sums = TRUE)

build_profile_pair_target_asset_v2(config, n_population = 100000L,
  repeats = 4L, seed = 87001L, h_target = NULL,
  development = FALSE, implementation_commit = "working_tree",
  control = list(), checkpoint = NULL)

build_profile_pair_direction_asset_v2(config, target_asset,
  n_population = 100000L, repeats = 4L, seed = 97001L,
  requested_coordinates = NULL, development = FALSE,
  implementation_commit = "working_tree", control = list(),
  checkpoint = NULL)

validate_profile_pair_target_asset_v2(asset, config = NULL,
  final = TRUE, expected_commit = NULL)

validate_profile_pair_direction_asset_v2(asset, config = NULL,
  final = TRUE, expected_commit = NULL)

profile_pair_target_scores_v2(dat, target_asset,
  h_analysis = target_asset$h_analysis,
  nuisance_control = list(grad_tol = 1e-7))
}
\arguments{
\item{y}{Finite numeric response vector, with one entry per observation.}
\item{X}{Finite numeric fixed-effect design matrix with unique nonempty column
names. Its rows correspond to \code{y}.}
\item{Z}{Finite numeric nuisance design matrix. A column of ones represents
a cluster intercept; an additional time column represents a cluster slope.}
\item{cluster_id}{Cluster identifier for every observation. Independent
clusters define the sample size used by the procedure.}
\item{tau}{Quantile level in \eqn{(0,1)}.}
\item{h_inf, h}{Positive inferential bandwidth. The two constituent profiles
use \eqn{h/2} and \eqn{h}. A prospective choice is \eqn{n^{-3/10}}, where
\eqn{n} is the number of clusters.}
\item{h_est}{Positive bandwidth for the convex selection fit. With
\code{NULL}, the rule is \eqn{\{\log(\max(p_P,2))/n\}^{1/4}}, where
\eqn{p_P} is the number of penalised columns.}
\item{lambda_gamma}{Fixed positive nuisance penalty. Supply a scalar,
a positive diagonal vector, or a symmetric positive-definite matrix
\eqn{\Lambda} of size \code{ncol(Z)}. The same matrix is used throughout
selection, refitting, and both inferential profiles.}
\item{mean_score_multiplier}{Positive prospective multiplier \eqn{c_b} of
the \code{h_est} penalty component. The default is one. Its theoretical
adequacy depends on the score domination conditions described below.}
\item{lambda_multiplier}{Positive prospective multiplier applied to the sum
of the stochastic and mean-score penalty components.}
\item{penalty_factor}{Nonnegative vector of length \code{ncol(X)}. Zero
marks an always-included coordinate. Positive entries multiply the cluster
score loadings; the default assigns one to every coordinate.}
\item{target_coordinates}{Prespecified coordinate names or integer indices
to retain in the refit. Declare every coordinate intended for
\code{infer_profile_pair_v2} when constructing the fit.}
\item{control}{For \code{fit_profile_pair_v2}, a named list containing
\code{fit_control}, \code{refit_control}, \code{nuisance_control}, or
\code{assembly}. The first two are optimiser controls; the nuisance controls
apply to every stage. For the population builders, controls include
\code{max_iter}, \code{max_root_iter}, \code{score_tol},
\code{nuisance_tol}, \code{identity_tol}, \code{eigen_min}, and
\code{condition_max}. Population acceptance limits remain enforced.}
\item{fit}{A converged object returned by \code{fit_profile_pair_v2}.}
\item{coordinates}{Unique names or integer indices among the coordinates
prespecified in \code{fit}.}
\item{mu_grid}{Requested positive finite Dantzig tolerances. The admissible
domain for coordinate inference is \eqn{0<\mu<1}. Requested values at least
one are retained with status \code{excluded_zero_direction}.}
\item{ci_level}{Pointwise confidence level in \eqn{(0,1)}.}
\item{solver_preference}{Requested solver names including \code{"CLARABEL"}.
The implemented paired path uses native CLARABEL and records that solver.}
\item{solver_opts}{Named arguments passed to the native CLARABEL solver.
Solver settings belong in its \code{control} list. Defaults include
\code{verbose = FALSE} and
\code{static_regularization_constant = 1e-6}. The returned
\code{solver_options} records the settings actually supplied.}
\item{df_correction}{Logical. Apply the prospective factor \eqn{n/(n-1)} to
the cluster variance when \code{TRUE}.}
\item{beta}{Finite coefficient vector at which the two complete profiles
are evaluated.}
\item{nuisance_control}{Named controls for the low-dimensional nuisance
optimiser, including \code{grad_tol} and \code{newton_maxit}.}
\item{identity_check}{Logical. Compute and compare the Schur and
ridge-curvature expressions for each constituent Hessian.}
\item{assembly}{Matrix assembly mode, \code{"sequential"} or
\code{"batch"}. \code{NULL} uses the configured mode.}
\item{need_fold_sums}{Logical. Retain deterministic cluster-fold Hessian
sums when \code{TRUE}.}
\item{config}{Named registered scenario list. It specifies the full
fixed-effect dimension \code{p}, sparsity \code{s}, quantile \code{tau},
nuisance dimension \code{q}, analysis cluster count \code{n_clusters},
cluster-size and design settings, error and cluster-effect distributions,
\code{fit_random_effects}, \code{lambda_gamma}, \code{h_multiplier},
\code{target_mode}, and \code{experiment_id}.}
\item{n_population}{Number of independently generated clusters per
population repeat. Final assets require at least 100000.}
\item{repeats}{Number of independent population repeats. Final assets
require at least four; development construction requires at least two.}
\item{seed}{Prospective seed used to derive deterministic role-specific
repeat seeds. Target and direction constructions use disjoint seeds.}
\item{h_target}{Positive target-approximation bandwidth smaller than the
analysis bandwidth. The default is the smaller of one quarter of the
analysis bandwidth and \eqn{n_{\mathrm{population}}^{-1/3}}.}
\item{development}{Logical. Permit development-scale construction when
\code{TRUE}. Its value is stored in the asset.}
\item{implementation_commit}{Character identifier of the implementation
used to construct the asset.}
\item{checkpoint}{Optional function called as \code{checkpoint(label, value)}
after a completed population step.}
\item{target_asset}{A paired profile target asset returned by
\code{build_profile_pair_target_asset_v2}.}
\item{requested_coordinates}{Unique coordinate names for population
precision directions. \code{NULL} requests all fixed-effect coordinates.}
\item{asset}{The paired target or direction asset to validate.}
\item{final}{Logical. Require final population scale and a final asset
when \code{TRUE}.}
\item{expected_commit}{Optional implementation identifier that must match
the stored asset identifier.}
\item{dat}{Clustered data list containing \code{y}, \code{X}, \code{Z},
\code{cluster_id}, and the variables needed by the asset's working nuisance
design. Column names must match the full target design.}
\item{h_analysis}{Positive bandwidth used for target-level score diagnostics.
Defaults to the analysis bandwidth recorded in \code{target_asset}.}
}
\details{
\subsection{Profile target and score}{
The inferential target minimises the population expectation of the complete
unsmoothed check-loss criterion after profiling a cluster-specific nuisance
vector with penalty \eqn{\gamma^T\Lambda\gamma/2}. At each supplied
coefficient vector, the nuisance optimisers are computed separately at the
two bandwidths. The inferential score and its exact Jacobian are
\deqn{G_{i,h}^{R}(\beta)=2g_{i,h/2}(\beta)-g_{i,h}(\beta),}
\deqn{J_{i,h}^{R}(\beta)=2H_{i,h/2}(\beta)-H_{i,h}(\beta).}
Their recorded score definition is
\code{profile_bandwidth_pair_2_half_minus_full}. The one-step estimator
and the cluster sandwich use this same score. Combining each cluster's
contributions before computing its covariance retains the covariance
between bandwidths.
}
\subsection{Initial fit and numerical checks}{
Let \eqn{p_P} count the penalised coordinates. The stochastic component is
\eqn{1.10\Phi^{-1}\{1-\alpha_n/(2p_P)\}/\sqrt n}, where
\eqn{\alpha_n=0.10/\log\{\max(n,3)\}}. The mean-score component is
\eqn{c_bh_{\mathrm{est}}}. Their sum, multiplied by
\code{lambda_multiplier} and the final coordinate loading, defines the
penalty. Loadings are computed at zero and at a preliminary fit; their
coordinatewise maximum is fixed during the final convex L1 optimisation.

Theoretical score domination must cover the centred stochastic score, the
population mean score of order \eqn{h_{\mathrm{est}}}, and numerical KKT
error. For example, with loading lower bound \eqn{w_{\min}>0}, a
sufficient scalar penalty level is
\deqn{\lambda_{\mathrm{est}}\geq
2\{C_g a_n+B h_{\mathrm{est}}+\eta_n\}/w_{\min},}
where \eqn{C_g a_n} bounds the centred score, \eqn{B h_{\mathrm{est}}}
bounds its population mean, and \eqn{\eta_n} bounds numerical KKT error.
The default \eqn{c_b=1} specifies the numerical procedure. Valid inference
also requires the stated model-class domination, positive loading bounds,
and effective profile curvature conditions.

The refit contains selected coordinates, always-included coordinates, and
declared inference targets. Its dimension must be at most
\eqn{\lfloor n/\log\{\max(n,3)\}\rfloor}. A zero-L1 convex profile fit
at \code{h_inf} supplies \code{beta_refit}. Both inferential profiles then
use that same full-design coefficient vector. Numerical failures terminate
with class \code{profile_pair_error}, carrying \code{failure_stage}, the
available fit record, and diagnostics.
}
\subsection{Admissible precision tolerance}{
For canonical \eqn{e_k}, the Dantzig constraint is
\eqn{\|J_{n,h}^{R}\omega-e_k\|_\infty\leq\mu}. At \eqn{\omega=0}
its left side is one. Every \eqn{\mu\geq1} therefore has the zero vector
as its unique L1 minimiser. Before cross-validation the procedure forms
the requested grid's intersection with \eqn{(0,1)} and records all excluded
values. An empty admissible grid terminates at \code{precision_solver}.

Cross-validation uses two deterministic cluster folds for \eqn{n<200}
and four otherwise. Each candidate is evaluated by the held-out inverse
defect using paired fold Jacobians. Candidates require at least half the
folds to return a defect. The smallest mean defect is selected, and numerical
ties select the smaller tolerance. The selected direction is solved on
the full paired Jacobian. Population-direction feasibility, precision
rates, and positive projected-score variance remain inferential conditions.
}
\subsection{Population assets}{
Target construction covers every registered fixed-effect column and records
independent repeat estimates at target, refinement, and analysis bandwidths.
The direction builder uses a separate population sample. It stores the
paired Hessian, its inverse, and the matched cluster-score covariance,
including the covariance between constituent bandwidths. Validators check
method identity, design coverage, convergence, and declared acceptance
conditions. Final validation additionally requires the final population
scale and \code{development=FALSE}.
}
}
\value{
\code{fit_profile_pair_v2} returns a \code{profile_pair_fit_v2} list containing
\code{beta_l1}, \code{beta_refit}, \code{beta}, paired \code{components},
\code{components_first_stage}, loading passes, \code{lambda_components},
selected and refit sets, numerical diagnostics, and method identifiers.

\code{infer_profile_pair_v2} returns a list with \code{table},
\code{directions}, \code{cv}, \code{fit}, and the requested and admissible
tolerance grids. Its table contains coordinate estimates, paired standard
errors, intervals, and solver diagnostics. \code{se_smoothed} is a
constituent-bandwidth diagnostic. The field \code{excluded_zero_direction}
records values outside the admissible tolerance domain.

\code{profile_bandwidth_pair_v2} returns the combined objective, score,
Jacobian, cluster scores, fold sums, and both constituent component lists.

The population builders return versioned lists with repeat records,
population estimates, Monte Carlo uncertainty, bandwidth identity, and
acceptance diagnostics. The target asset includes \code{beta_star_mc};
the direction asset includes \code{H_population}, \code{Omega_population},
\code{Sigma_population}, and named \code{directions}.

The validators invisibly return
\code{list(valid = TRUE, problems = character())} after successful checks
and terminate on an incompatible record. \code{profile_pair_target_scores_v2}
returns the paired cluster scores, mean score, Jacobian, bandwidths, and
geometry diagnostics evaluated at the target asset's coefficient vector.
}
\examples{
set.seed(20261022)
n <- 20L
m <- 4L
p <- 5L
cluster_id <- rep(seq_len(n), each = m)
X <- matrix(rnorm(n * m * p), n * m, p,
            dimnames = list(NULL, paste0("x", seq_len(p))))
Z <- matrix(1, n * m, 1L)
beta <- c(0.75, -0.5, 0, 0, 0)
b <- rnorm(n)
y <- as.numeric(X \%*\% beta) + rep(b, each = m) + rnorm(n * m)
fit <- fit_profile_pair_v2(
  y, X, Z, cluster_id, tau = 0.5, h_inf = n^(-3 / 10),
  target_coordinates = c("x1", "x4")
)
inference <- infer_profile_pair_v2(
  fit, c("x1", "x4"), mu_grid = c(0.05, 0.15, 0.4)
)
inference$table[, c("coordinate", "beta_tilde", "se")]
}
\keyword{regression}
