#include <RcppEigen.h>
#include <Highs.h>
#include <algorithm>
#include <cmath>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp17)]]

namespace {
void soy_highs_status(HighsStatus status) {
  if (status == HighsStatus::kError) Rcpp::stop("HiGHS rejected a pooled quantile operation.");
}
}

// [[Rcpp::export]]
Rcpp::List soy_pooled_lp_cpp(const Eigen::Map<Eigen::MatrixXd> X,
                            const Eigen::Map<Eigen::VectorXd> y,
                            const Eigen::Map<Eigen::VectorXd> factor,
                            double tau, int nlambda, double lambda_min,
                            Rcpp::NumericVector requested, bool progress) {
  const int N = X.rows(), p = X.cols();
  if (y.size() != N || factor.size() != p || !X.allFinite() || !y.allFinite() ||
      tau <= 0 || tau >= 1 || nlambda < 2 || lambda_min <= 0 || lambda_min > 1)
    Rcpp::stop("Invalid pooled quantile path input.");
  Highs solver;
  soy_highs_status(solver.setOptionValue("output_flag", false));
  soy_highs_status(solver.setOptionValue("threads", 1));
  soy_highs_status(solver.setOptionValue("solver", "simplex"));
  soy_highs_status(solver.setOptionValue("presolve", "off"));
  soy_highs_status(solver.setOptionValue("simplex_strategy", 1));
  soy_highs_status(solver.setOptionValue("simplex_scale_strategy", 2));
  soy_highs_status(solver.setOptionValue("primal_feasibility_tolerance", 1e-9));
  soy_highs_status(solver.setOptionValue("dual_feasibility_tolerance", 1e-9));
  HighsLp model;
  model.num_col_ = N;
  model.num_row_ = 0;
  model.sense_ = ObjSense::kMaximize;
  model.col_cost_.assign(y.data(), y.data() + N);
  model.col_lower_.assign(N, tau - 1);
  model.col_upper_.assign(N, tau);
  model.a_matrix_.format_ = MatrixFormat::kColwise;
  model.a_matrix_.start_.assign(N + 1, 0);
  soy_highs_status(solver.passModel(model));
  std::vector<int> selected;
  std::vector<bool> present(p, false);
  auto add = [&](const std::vector<int>& columns, double lambda) {
    std::vector<HighsInt> starts(1, 0), indices;
    std::vector<double> values, lower, upper;
    for (int k : columns) {
      for (int i = 0; i < N; ++i) if (X(i, k) != 0) {
        indices.push_back(i);
        values.push_back(X(i, k));
      }
      starts.push_back(indices.size());
      lower.push_back(-N * lambda * factor[k]);
      upper.push_back(N * lambda * factor[k]);
      selected.push_back(k);
      present[k] = true;
    }
    soy_highs_status(solver.addRows(columns.size(), lower.data(), upper.data(), values.size(),
      starts.data(), indices.data(), values.data()));
  };
  std::vector<int> fixed;
  for (int k = 0; k < p; ++k) {
    if (factor[k] == 0) fixed.push_back(k);
    else if (factor[k] != 1) Rcpp::stop("Pooled penalty factors must be zero or one.");
  }
  if (fixed.empty() || fixed.size() == static_cast<size_t>(p)) Rcpp::stop("The pooled design requires both adjustment and penalised columns.");
  add(fixed, 0);
  auto solve = [&]() {
    Rcpp::checkUserInterrupt();
    soy_highs_status(solver.run());
    if (solver.getModelStatus() != HighsModelStatus::kOptimal)
      Rcpp::stop("A pooled quantile LP did not reach its optimum.");
  };
  solve();
  Eigen::VectorXd dual = Eigen::Map<const Eigen::VectorXd>(solver.getSolution().col_value.data(), N);
  Eigen::VectorXd gradient = X.transpose() * dual;
  std::vector<double> lambdas;
  if (requested.size()) {
    lambdas.assign(requested.begin(), requested.end());
  } else {
    double maximum = 0;
    for (int k = 0; k < p; ++k) if (factor[k]) maximum = std::max(maximum, std::abs(gradient[k]) / N);
    if (!std::isfinite(maximum) || maximum <= 0) Rcpp::stop("The pooled entry penalty is not positive.");
    for (int j = 0; j < nlambda; ++j)
      lambdas.push_back(maximum * std::exp(std::log(lambda_min) * j / (nlambda - 1)));
  }
  const int count = lambdas.size();
  Eigen::MatrixXd coefficients = Eigen::MatrixXd::Zero(p, count);
  Rcpp::NumericVector gaps(count), violations(count);
  Rcpp::IntegerVector constraints(count), nonzero(count), generations(count), iterations(count);
  for (int point = 0; point < count; ++point) {
    const double lambda = lambdas[point];
    if (!std::isfinite(lambda) || lambda < 0 || (point && lambda > lambdas[point - 1]))
      Rcpp::stop("Invalid decreasing pooled penalty grid.");
    std::vector<double> lower, upper;
    for (int k : selected) {
      lower.push_back(-N * lambda * factor[k]);
      upper.push_back(N * lambda * factor[k]);
    }
    soy_highs_status(solver.changeRowsBounds(0, selected.size() - 1, lower.data(), upper.data()));
    int generation = 0, simplex_iterations = 0;
    while (true) {
      if (++generation > p) Rcpp::stop("Pooled constraint generation exhausted its full coordinate budget.");
      solve();
      simplex_iterations += solver.getInfo().simplex_iteration_count;
      dual = Eigen::Map<const Eigen::VectorXd>(solver.getSolution().col_value.data(), N);
      gradient.noalias() = X.transpose() * dual;
      std::vector<int> additions;
      for (int k = 0; k < p; ++k)
        if (!present[k] && (std::abs(gradient[k]) - N * lambda * factor[k]) / N > 1e-9)
          additions.push_back(k);
      if (additions.empty()) break;
      std::stable_sort(additions.begin(), additions.end(), [&](int a, int b) {
        return std::abs(gradient[a]) - N * lambda * factor[a] > std::abs(gradient[b]) - N * lambda * factor[b];
      });
      if (additions.size() > 50) additions.resize(50);
      add(additions, lambda);
    }
    const auto& solution = solver.getSolution();
    Eigen::VectorXd beta = Eigen::VectorXd::Zero(p);
    for (size_t r = 0; r < selected.size(); ++r) beta[selected[r]] = solution.row_dual[r];
    const Eigen::VectorXd residual = y - X * beta;
    double primal = 0, violation = 0;
    for (int i = 0; i < N; ++i) primal += residual[i] * (tau - (residual[i] < 0));
    for (int k = 0; k < p; ++k) {
      primal += N * lambda * factor[k] * std::abs(beta[k]);
      violation = std::max(violation, (std::abs(gradient[k]) - N * lambda * factor[k]) / N);
      if (std::abs(beta[k]) > 1e-10) ++nonzero[point];
    }
    const double gap = (primal - y.dot(dual)) / N;
    if (!std::isfinite(gap) || std::abs(gap) > 1e-7 || violation > 1e-8 ||
        dual.minCoeff() < tau - 1 - 1e-8 || dual.maxCoeff() > tau + 1e-8)
      Rcpp::stop("The complete pooled primal-dual certificate failed.");
    coefficients.col(point) = beta;
    gaps[point] = gap;
    violations[point] = violation;
    constraints[point] = selected.size();
    generations[point] = generation;
    iterations[point] = simplex_iterations;
    if (progress) Rcpp::Rcout << "pooled penalty " << point + 1 << " of " << count
      << " constraints " << selected.size() << " gap " << gap << std::endl;
    if (point + 1 < count) {
      const HighsBasis basis = solver.getBasis();
      std::vector<HighsInt> removed;
      std::vector<int> retained;
      for (size_t r = 0; r < selected.size(); ++r) {
        const int k = selected[r];
        if (factor[k] && basis.row_status[r] == HighsBasisStatus::kBasic) {
          removed.push_back(r);
          present[k] = false;
        } else retained.push_back(k);
      }
      if (!removed.empty()) {
        soy_highs_status(solver.deleteRows(removed.size(), removed.data()));
        selected.swap(retained);
      }
    }
  }
  return Rcpp::List::create(Rcpp::_["beta"] = coefficients, Rcpp::_["lambda"] = lambdas,
    Rcpp::_["primal_dual_gap"] = gaps, Rcpp::_["dual_violation"] = violations,
    Rcpp::_["active_constraints"] = constraints, Rcpp::_["selected_variables"] = nonzero,
    Rcpp::_["generations"] = generations, Rcpp::_["iterations"] = iterations,
    Rcpp::_["highs_version"] = solver.version());
}
