#include <RcppArmadillo.h>
#include <R_ext/Rdynload.h>

namespace {
double option_number(const Rcpp::List& control, const char* name, double value) {
  return control.containsElementNamed(name) ? Rcpp::as<double>(control[name]) : value;
}

double smooth_loss(double residual, double tau, double bandwidth) {
  const double cut = std::max(-1.0, std::min(1.0, residual / bandwidth));
  const double distribution = cut <= -1.0 ? 0.0 :
    (cut >= 1.0 ? 1.0 : 0.5 + 0.75 * (cut - cut * cut * cut / 3.0));
  const double primitive = 0.375 * cut * cut - 0.1875 * std::pow(cut, 4.0);
  const double first = primitive - 0.1875;
  return tau * (residual * distribution - bandwidth * first) +
    (tau - 1.0) * (residual * (1.0 - distribution) + bandwidth * first);
}

void evaluate(const arma::vec& outcome, const arma::mat& design,
              const arma::vec& gamma, const arma::mat& penalty,
              double tau, double bandwidth, double& objective,
              arma::vec& gradient, arma::mat& hessian,
              arma::vec& residual, arma::vec& score, arma::vec& curvature) {
  residual = outcome - design * gamma;
  score.set_size(residual.n_elem);
  curvature.set_size(residual.n_elem);
  double loss = 0.0;
  for (arma::uword j = 0; j < residual.n_elem; ++j) {
    const double u = residual[j] / bandwidth;
    const double v = -u;
    const double distribution = v <= -1.0 ? 0.0 :
      (v >= 1.0 ? 1.0 : 0.5 + 0.75 * (v - v * v * v / 3.0));
    score[j] = tau - distribution;
    curvature[j] = std::abs(u) <= 1.0 ? 0.75 * (1.0 - u * u) / bandwidth : 0.0;
    loss += smooth_loss(residual[j], tau, bandwidth);
  }
  const double size = static_cast<double>(residual.n_elem);
  objective = loss / size + 0.5 * arma::dot(gamma, penalty * gamma);
  gradient = -(design.t() * score) / size + penalty * gamma;
  arma::mat weighted = design;
  weighted.each_col() %= curvature / size;
  hessian = design.t() * weighted + penalty;
}
}

extern "C" SEXP jpdmeqr_cluster_profiles(SEXP outcome_sexp, SEXP design_sexp,
                                        SEXP rows_sexp, SEXP pointers_sexp,
                                        SEXP initial_sexp, SEXP tau_sexp,
                                        SEXP bandwidth_sexp, SEXP penalty_sexp,
                                        SEXP control_sexp) {
  BEGIN_RCPP
  const arma::vec outcome = Rcpp::as<arma::vec>(outcome_sexp);
  const arma::mat design = Rcpp::as<arma::mat>(design_sexp);
  const Rcpp::IntegerVector rows(rows_sexp), pointers(pointers_sexp);
  const arma::mat initial = Rcpp::as<arma::mat>(initial_sexp);
  const arma::mat penalty = Rcpp::as<arma::mat>(penalty_sexp);
  const Rcpp::List control(control_sexp);
  const double tau = Rcpp::as<double>(tau_sexp);
  const double bandwidth = Rcpp::as<double>(bandwidth_sexp);
  const double tolerance = option_number(control, "grad_tol", 1e-8);
  const int max_iterations = static_cast<int>(option_number(control, "newton_maxit", 100));
  const int max_backtracks = static_cast<int>(option_number(control, "max_backtrack", 50));
  const double armijo = option_number(control, "armijo", 1e-4);
  const double reduction = option_number(control, "backtrack_factor", 0.5);
  const double minimum_step = option_number(control, "min_step", 1e-14);
  const int clusters = pointers.size() - 1;
  if (clusters < 1 || rows.size() != static_cast<int>(outcome.n_elem) ||
      design.n_rows != outcome.n_elem || design.n_cols != penalty.n_rows ||
      penalty.n_rows != penalty.n_cols || initial.n_rows != static_cast<arma::uword>(clusters) ||
      initial.n_cols != design.n_cols || pointers[0] != 0 ||
      pointers[clusters] != rows.size() || !(tau > 0 && tau < 1) || !(bandwidth > 0) ||
      !(tolerance > 0) || max_iterations < 1 || max_backtracks < 0 ||
      !(reduction > 0 && reduction < 1) || !(armijo > 0 && armijo < 1) ||
      !(minimum_step > 0) || !outcome.is_finite() || !design.is_finite() ||
      !initial.is_finite() || !penalty.is_finite()) {
    Rcpp::stop("Invalid compiled cluster-profile input.");
  }
  arma::mat gamma_matrix(clusters, design.n_cols), gradients(clusters, design.n_cols);
  arma::vec objectives(clusters), gradient_max(clusters), residual_all(outcome.n_elem),
    score_all(outcome.n_elem), curvature_all(outcome.n_elem);
  Rcpp::IntegerVector iterations(clusters), backtracks(clusters), evaluations(clusters);
  for (int i = 0; i < clusters; ++i) {
    if (i % 128 == 0) Rcpp::checkUserInterrupt();
    const int begin = pointers[i], end = pointers[i + 1];
    if (begin < 0 || end <= begin || end > rows.size()) Rcpp::stop("Invalid cluster pointers.");
    arma::uvec index(end - begin);
    for (int j = begin; j < end; ++j) {
      if (rows[j] < 0 || rows[j] >= static_cast<int>(outcome.n_elem)) {
        Rcpp::stop("Invalid cluster row index.");
      }
      index[j - begin] = rows[j];
    }
    const arma::vec cluster_outcome = outcome.elem(index);
    const arma::mat cluster_design = design.rows(index);
    arma::vec gamma = initial.row(i).t(), gradient, residual, score, curvature;
    arma::mat hessian;
    double objective;
    evaluate(cluster_outcome, cluster_design, gamma, penalty, tau, bandwidth,
             objective, gradient, hessian, residual, score, curvature);
    evaluations[i] = 1;
    for (int iteration = 0; iteration < max_iterations; ++iteration) {
      if (arma::abs(gradient).max() <= tolerance) break;
      arma::vec direction;
      const bool solved = arma::solve(direction, hessian, -gradient,
        arma::solve_opts::likely_sympd + arma::solve_opts::no_approx);
      const double slope = solved ? arma::dot(gradient, direction) : 0.0;
      if (!solved || !direction.is_finite() || !(slope < 0)) {
        Rcpp::stop("Compiled nuisance Newton direction failed at cluster %d.", i + 1);
      }
      double step = 1.0;
      bool accepted = false;
      arma::vec candidate, candidate_gradient, candidate_residual, candidate_score, candidate_curvature;
      arma::mat candidate_hessian;
      double candidate_objective = NA_REAL;
      for (int attempt = 0; attempt <= max_backtracks; ++attempt) {
        candidate = gamma + step * direction;
        evaluate(cluster_outcome, cluster_design, candidate, penalty, tau, bandwidth,
                 candidate_objective, candidate_gradient, candidate_hessian,
                 candidate_residual, candidate_score, candidate_curvature);
        evaluations[i] += 1;
        const double rounding = 16.0 * std::numeric_limits<double>::epsilon() *
          std::max(1.0, std::abs(objective));
        if (std::isfinite(candidate_objective) &&
            candidate_objective <= objective + armijo * step * slope + rounding) {
          accepted = true;
          backtracks[i] += attempt;
          break;
        }
        step *= reduction;
        if (step < minimum_step) break;
      }
      if (!accepted) Rcpp::stop("Compiled nuisance line search failed at cluster %d.", i + 1);
      gamma = candidate;
      gradient = candidate_gradient;
      hessian = candidate_hessian;
      residual = candidate_residual;
      score = candidate_score;
      curvature = candidate_curvature;
      objective = candidate_objective;
      iterations[i] = iteration + 1;
    }
    const double maximum = arma::abs(gradient).max();
    if (!std::isfinite(objective) || !gradient.is_finite() || maximum > tolerance) {
      Rcpp::stop("Compiled nuisance tolerance failed at cluster %d: %.17g.", i + 1, maximum);
    }
    gamma_matrix.row(i) = gamma.t();
    gradients.row(i) = gradient.t();
    objectives[i] = objective;
    gradient_max[i] = maximum;
    residual_all.elem(index) = residual;
    score_all.elem(index) = score;
    curvature_all.elem(index) = curvature;
  }
  return Rcpp::List::create(
    Rcpp::_["gamma"] = gamma_matrix, Rcpp::_["gradient"] = gradients,
    Rcpp::_["objective"] = objectives, Rcpp::_["gradient_max"] = gradient_max,
    Rcpp::_["residual"] = residual_all, Rcpp::_["psi"] = score_all,
    Rcpp::_["phi"] = curvature_all, Rcpp::_["iterations"] = iterations,
    Rcpp::_["backtracks"] = backtracks, Rcpp::_["evaluations"] = evaluations);
  END_RCPP
}

extern "C" SEXP jpDMEQR_population_crossprod(SEXP x, SEXP y);

static const R_CallMethodDef CallEntries[] = {
  {"jpdmeqr_cluster_profiles", reinterpret_cast<DL_FUNC>(&jpdmeqr_cluster_profiles), 9},
  {"jpDMEQR_population_crossprod", reinterpret_cast<DL_FUNC>(&jpDMEQR_population_crossprod), 2},
  {NULL, NULL, 0}
};

extern "C" void R_init_jpDMEQR(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
