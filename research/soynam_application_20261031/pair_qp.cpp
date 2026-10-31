#include <RcppEigen.h>
#include <vector>
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp17)]]

using Eigen::MatrixXd;
using Eigen::VectorXd;

// [[Rcpp::export]]
Rcpp::List soy_pair_qp_cpp(
    const Eigen::Map<MatrixXd> X, const Eigen::Map<MatrixXd> Z,
    const Eigen::Map<VectorXd> outcome, Rcpp::IntegerVector cluster,
    const Eigen::Map<VectorXd> inverse_size, const Eigen::Map<MatrixXd> Lambda,
    const Eigen::Map<VectorXd> coordinate_penalty, double quantile,
    VectorXd theta, MatrixXd gamma, VectorXd dual,
    int max_iterations = 100000, double gap_tolerance = 1e-8) {
  const int N = X.rows(), d = X.cols(), n = gamma.rows();
  if (Z.cols() != 1 || Lambda.rows() != 1 || Lambda.cols() != 1 || N != 2 * n ||
      !(Lambda(0, 0) > 0) || coordinate_penalty.cwiseAbs().maxCoeff() != 0 ||
      (Z.array() - 1).abs().maxCoeff() > 1e-12)
    Rcpp::stop("The selected solver requires paired intercept clusters and zero L1 penalty.");
  const double lambda = Lambda(0, 0), rho = 1.0, relaxation = 1.6;
  std::vector<int> first(n, -1), second(n, -1), previous_pattern(N, 9);
  MatrixXd C(d, n);
  for (int j = 0; j < N; ++j) {
    const int i = cluster[j];
    if (i < 0 || i >= n || std::abs(inverse_size[j] - 0.5) > 1e-12)
      Rcpp::stop("Invalid paired cluster weights.");
    if (first[i] < 0) first[i] = j;
    else if (second[i] < 0) second[i] = j;
    else Rcpp::stop("A cluster has too many observations.");
  }
  for (int i = 0; i < n; ++i) {
    if (second[i] < 0) Rcpp::stop("A paired observation is absent.");
    C.col(i) = 0.5 * (X.row(first[i]) + X.row(second[i])).transpose();
  }
  MatrixXd H = MatrixXd::Identity(d, d) + (0.5 * rho) * X.transpose() * X;
  H.noalias() -= (rho * rho / (lambda + rho)) * C * C.transpose();
  Eigen::LLT<MatrixXd> factor(H);
  if (factor.info() != Eigen::Success) Rcpp::stop("The paired quadratic system is not positive definite.");
  VectorXd effect = gamma.col(0);
  VectorXd lower = (quantile - 1) * inverse_size, upper = quantile * inverse_size;
  dual = dual.cwiseMax(lower).cwiseMin(upper);
  VectorXd u = -dual.array() / (rho * inverse_size.array());
  VectorXd prediction = X * theta;
  for (int j = 0; j < N; ++j) prediction[j] += effect[cluster[j]];
  VectorXd residual = outcome - prediction;
  double primal = R_PosInf, dual_value = R_NegInf, gap = R_PosInf;
  int face_checks = 0;
  auto evaluate = [&](const VectorXd& coefficient, const VectorXd& nuisance, const VectorXd& score,
                      double& p, double& v) {
    VectorXd r = outcome - X * coefficient;
    double loss = 0;
    for (int j = 0; j < N; ++j) {
      r[j] -= nuisance[cluster[j]];
      loss += inverse_size[j] * r[j] * (quantile - (r[j] < 0));
    }
    p = (loss + 0.5 * coefficient.squaredNorm() + 0.5 * lambda * nuisance.squaredNorm()) / n;
    VectorXd conjugate = X.transpose() * score;
    double gamma_conjugate = 0;
    for (int i = 0; i < n; ++i) {
      const double value = score[first[i]] + score[second[i]];
      gamma_conjugate += value * value / lambda;
    }
    v = (outcome.dot(score) - 0.5 * conjugate.squaredNorm() - 0.5 * gamma_conjugate) / n;
    return p - v;
  };
  auto finish = [&](int iteration, bool polished) {
    gamma = effect;
    VectorXd beta_gradient = theta - X.transpose() * dual;
    double gamma_error = 0;
    for (int i = 0; i < n; ++i)
      gamma_error = std::max(gamma_error, std::abs(lambda * effect[i] - dual[first[i]] - dual[second[i]]));
    return Rcpp::List::create(Rcpp::_["theta"] = theta, Rcpp::_["gamma"] = gamma,
      Rcpp::_["dual"] = dual, Rcpp::_["primal"] = primal, Rcpp::_["dual_value"] = dual_value,
      Rcpp::_["duality_gap"] = gap, Rcpp::_["beta_kkt"] = beta_gradient.cwiseAbs().maxCoeff() / n,
      Rcpp::_["gamma_kkt"] = gamma_error, Rcpp::_["iterations"] = iteration,
      Rcpp::_["face_checks"] = face_checks, Rcpp::_["polished"] = polished,
      Rcpp::_["solver"] = "paired_structured_qp", Rcpp::_["solver_status"] = "certified");
  };
  for (int iteration = 1; iteration <= max_iterations; ++iteration) {
    Rcpp::checkUserInterrupt();
    VectorXd v = outcome - residual - u, means(n);
    for (int i = 0; i < n; ++i) means[i] = 0.5 * (v[first[i]] + v[second[i]]);
    VectorXd rhs = 0.5 * rho * X.transpose() * v - (rho * rho / (lambda + rho)) * C * means;
    theta = factor.solve(rhs);
    effect = (rho / (lambda + rho)) * (means - C.transpose() * theta);
    prediction.noalias() = X * theta;
    for (int j = 0; j < N; ++j) prediction[j] += effect[cluster[j]];
    VectorXd relaxed = relaxation * prediction + (1 - relaxation) * (outcome - residual);
    VectorXd argument = outcome - relaxed - u;
    for (int j = 0; j < N; ++j) {
      if (argument[j] > quantile / rho) residual[j] = argument[j] - quantile / rho;
      else if (argument[j] < (quantile - 1) / rho) residual[j] = argument[j] - (quantile - 1) / rho;
      else residual[j] = 0;
    }
    u += relaxed + residual - outcome;
    dual = (-rho * inverse_size.array() * u.array()).matrix();
    dual = dual.cwiseMax(lower).cwiseMin(upper);
    gap = evaluate(theta, effect, dual, primal, dual_value);
    if (!std::isfinite(gap) || gap < -1e-10 * std::max(1.0, std::abs(primal)))
      Rcpp::stop("Invalid paired quadratic certificate.");
    if (gap <= gap_tolerance * std::max(1.0, std::abs(primal))) return finish(iteration, false);
    if (iteration != 1 && iteration % 5 != 0) continue;
    std::vector<int> pattern(N), active_rows, double_clusters;
    bool changed = false;
    for (int j = 0; j < N; ++j) {
      pattern[j] = residual[j] > 0 ? 1 : (residual[j] < 0 ? -1 : 0);
      changed = changed || pattern[j] != previous_pattern[j];
    }
    if (!changed) continue;
    previous_pattern = pattern;
    VectorXd right = VectorXd::Zero(d);
    for (int i = 0; i < n; ++i) {
      const int a = first[i], b = second[i];
      if (pattern[a] == 0 && pattern[b] == 0) {
        active_rows.push_back(a);
        double_clusters.push_back(i);
        right.noalias() += lambda * outcome[a] * X.row(a).transpose();
      } else if (pattern[a] == 0 || pattern[b] == 0) {
        const int active = pattern[a] == 0 ? a : b;
        const int other = active == a ? b : a;
        active_rows.push_back(active);
        const double score = pattern[other] > 0 ? upper[other] : lower[other];
        right.noalias() += lambda * outcome[active] * X.row(active).transpose()
          + score * (X.row(other) - X.row(active)).transpose();
      } else {
        right.noalias() += (pattern[a] > 0 ? upper[a] : lower[a]) * X.row(a).transpose()
          + (pattern[b] > 0 ? upper[b] : lower[b]) * X.row(b).transpose();
      }
    }
    if (static_cast<int>(double_clusters.size()) > d) continue;
    ++face_checks;
    MatrixXd active_design(active_rows.size(), d);
    for (unsigned int j = 0; j < active_rows.size(); ++j) active_design.row(j) = X.row(active_rows[j]);
    MatrixXd face_H = MatrixXd::Identity(d, d) + lambda * active_design.transpose() * active_design;
    Eigen::LLT<MatrixXd> face_factor(face_H);
    if (face_factor.info() != Eigen::Success) Rcpp::stop("Invalid active-face quadratic system.");
    VectorXd candidate = face_factor.solve(right), multipliers(double_clusters.size());
    if (!double_clusters.empty()) {
      MatrixXd constraints(double_clusters.size(), d);
      VectorXd values(double_clusters.size());
      for (unsigned int j = 0; j < double_clusters.size(); ++j) {
        const int i = double_clusters[j];
        constraints.row(j) = X.row(first[i]) - X.row(second[i]);
        values[j] = outcome[first[i]] - outcome[second[i]];
      }
      MatrixXd inverse_constraints = face_factor.solve(constraints.transpose());
      MatrixXd schur = constraints * inverse_constraints;
      Eigen::ColPivHouseholderQR<MatrixXd> decomposition(schur);
      decomposition.setThreshold(1e-12);
      if (decomposition.rank() != schur.rows()) continue;
      multipliers = decomposition.solve(constraints * candidate - values);
      candidate.noalias() -= inverse_constraints * multipliers;
    }
    VectorXd candidate_effect(n), candidate_dual(N);
    unsigned int double_index = 0;
    for (int i = 0; i < n; ++i) {
      const int a = first[i], b = second[i];
      if (pattern[a] == 0 && pattern[b] == 0) {
        candidate_effect[i] = outcome[a] - X.row(a).dot(candidate);
        candidate_dual[b] = multipliers[double_index++];
        candidate_dual[a] = lambda * candidate_effect[i] - candidate_dual[b];
      } else if (pattern[a] == 0 || pattern[b] == 0) {
        const int active = pattern[a] == 0 ? a : b;
        const int other = active == a ? b : a;
        candidate_effect[i] = outcome[active] - X.row(active).dot(candidate);
        candidate_dual[other] = pattern[other] > 0 ? upper[other] : lower[other];
        candidate_dual[active] = lambda * candidate_effect[i] - candidate_dual[other];
      } else {
        candidate_dual[a] = pattern[a] > 0 ? upper[a] : lower[a];
        candidate_dual[b] = pattern[b] > 0 ? upper[b] : lower[b];
        candidate_effect[i] = (candidate_dual[a] + candidate_dual[b]) / lambda;
      }
    }
    candidate_dual = candidate_dual.cwiseMax(lower).cwiseMin(upper);
    double candidate_primal, candidate_value;
    const double candidate_gap = evaluate(candidate, candidate_effect, candidate_dual, candidate_primal, candidate_value);
    if (std::isfinite(candidate_gap) && candidate_gap >= -1e-10 * std::max(1.0, std::abs(candidate_primal)) &&
        candidate_gap <= gap_tolerance * std::max(1.0, std::abs(candidate_primal))) {
      theta.swap(candidate);
      effect.swap(candidate_effect);
      dual.swap(candidate_dual);
      primal = candidate_primal;
      dual_value = candidate_value;
      gap = candidate_gap;
      return finish(iteration, true);
    }
  }
  Rcpp::stop("The paired quadratic solver exhausted its certified iteration budget.");
}
