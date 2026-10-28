#include <RcppEigen.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>
// [[Rcpp::depends(RcppEigen)]]

using Eigen::MatrixXd;
using Eigen::VectorXd;

double threshold_value(double value, double threshold) {
  return value > threshold ? value - threshold :
    (value < -threshold ? value + threshold : 0.0);
}

double scalar_profile_root(const VectorXd& response, const VectorXd& design,
                           double lambda, double quantile) {
  const int m = response.size();
  std::vector<std::pair<double, double> > knots;
  double total = 0.0;
  for (int j = 0; j < m; ++j) {
    if (design[j] == 0.0) continue;
    const double weight = std::abs(design[j]) / m;
    knots.emplace_back(response[j] / design[j], weight);
    total += weight * (design[j] > 0 ? quantile : 1 - quantile);
  }
  std::sort(knots.begin(), knots.end());
  double cumulative = 0.0;
  for (std::size_t k = 0; k < knots.size();) {
    const double location = knots[k].first;
    const double candidate = (total - cumulative) / lambda;
    if (candidate <= location) return(candidate);
    double jump = 0.0;
    while (k < knots.size() && knots[k].first == location) jump += knots[k++].second;
    if (lambda * location + cumulative + jump - total >= 0) return(location);
    cumulative += jump;
  }
  return((total - cumulative) / lambda);
}

// [[Rcpp::export]]
Rcpp::List exact_cluster_profile_cpp(
    const Eigen::Map<VectorXd> outcome, const Eigen::Map<MatrixXd> Z,
    Rcpp::IntegerVector row_order, Rcpp::IntegerVector pointers,
    const Eigen::Map<MatrixXd> Lambda, double quantile,
    double tolerance = 1e-10, int max_iterations = 30000) {
  const int n = pointers.size() - 1, q = Z.cols(), N = outcome.size();
  if (n < 1 || q < 1 || Z.rows() != N || row_order.size() != N ||
      Lambda.rows() != q || Lambda.cols() != q || !(quantile > 0 && quantile < 1)) {
    Rcpp::stop("Invalid exact cluster input.");
  }
  Eigen::LLT<MatrixXd> factor(Lambda);
  if (factor.info() != Eigen::Success) Rcpp::stop("Lambda must be positive definite.");
  MatrixXd gamma = MatrixXd::Zero(n, q);
  VectorXd dual = VectorXd::Zero(N), losses(n), gaps(n);
  Rcpp::IntegerVector iterations(n);
  for (int i = 0; i < n; ++i) {
    const int begin = pointers[i], end = pointers[i + 1], m = end - begin;
    if (m < 1) Rcpp::stop("Every cluster must be nonempty.");
    MatrixXd Zi(m, q);
    VectorXd response(m);
    for (int j = 0; j < m; ++j) {
      const int row = row_order[begin + j];
      if (row < 0 || row >= N) Rcpp::stop("Invalid exact cluster row.");
      Zi.row(j) = Z.row(row);
      response[j] = outcome[row];
    }
    if (q == 1) {
      const double lambda = Lambda(0, 0);
      std::vector<std::pair<double, double> > knots;
      double total = 0.0;
      for (int j = 0; j < m; ++j) {
        const double z = Zi(j, 0);
        if (z == 0.0) continue;
        const double weight = std::abs(z) / m;
        knots.emplace_back(response[j] / z, weight);
        total += weight * (z > 0 ? quantile : 1.0 - quantile);
      }
      std::sort(knots.begin(), knots.end());
      double cumulative = 0.0, optimum = total / lambda;
      int events = 0;
      for (std::size_t k = 0; k < knots.size();) {
        const double location = knots[k].first;
        const double candidate = (total - cumulative) / lambda;
        ++events;
        if (candidate <= location) {
          optimum = candidate;
          break;
        }
        double jump = 0.0;
        while (k < knots.size() && knots[k].first == location) {
          jump += knots[k].second;
          ++k;
        }
        if (lambda * location + cumulative + jump - total >= 0.0) {
          optimum = location;
          break;
        }
        cumulative += jump;
        optimum = (total - cumulative) / lambda;
      }
      VectorXd u(m), residual = response - Zi.col(0) * optimum;
      std::vector<int> active;
      double stationarity = 0.0, active_lower_sum = 0.0;
      for (int j = 0; j < m; ++j) {
        const double roundoff = 64 * std::numeric_limits<double>::epsilon() *
          std::max(1.0, std::max(std::abs(response[j]), std::abs(Zi(j, 0) * optimum)));
        if (std::abs(residual[j]) <= roundoff && Zi(j, 0) != 0.0) {
          active.push_back(j);
          active_lower_sum += std::min(Zi(j, 0) * (quantile - 1) / m, Zi(j, 0) * quantile / m);
          u[j] = 0.0;
        } else {
          u[j] = residual[j] < 0 ? (quantile - 1) / m : quantile / m;
          stationarity += Zi(j, 0) * u[j];
        }
      }
      double remaining = lambda * optimum - stationarity - active_lower_sum;
      for (int j : active) {
        const double low = std::min(Zi(j, 0) * (quantile - 1) / m, Zi(j, 0) * quantile / m);
        const double high = std::max(Zi(j, 0) * (quantile - 1) / m, Zi(j, 0) * quantile / m);
        const double increment = std::max(0.0, std::min(high - low, remaining));
        u[j] = (low + increment) / Zi(j, 0);
        remaining -= increment;
      }
      double check = 0.0;
      for (int j = 0; j < m; ++j) check += residual[j] * (quantile - (residual[j] < 0 ? 1.0 : 0.0)) / m;
      const double primal = check + 0.5 * lambda * optimum * optimum;
      const double dual_nuisance = Zi.col(0).dot(u);
      const double dual_value = response.dot(u) - 0.5 * dual_nuisance * dual_nuisance / lambda;
      const double gap = primal - dual_value;
      if (!std::isfinite(gap) || gap < -1e-9 || gap > tolerance * std::max(1.0, std::abs(primal))) {
        Rcpp::stop("The scalar cluster profile failed its exact gap certificate.");
      }
      gamma(i, 0) = optimum;
      for (int j = 0; j < m; ++j) dual[row_order[begin + j]] = u[j];
      losses[i] = primal;
      gaps[i] = gap;
      iterations[i] = events;
      continue;
    }
    if (q == 2) {
      const MatrixXd transform = factor.matrixU().solve(MatrixXd::Identity(2, 2));
      const MatrixXd design = Zi * transform;
      const double low = (quantile - 1) / m, high = quantile / m;
      VectorXd solution_theta(2), solution_dual(m);
      double solution_loss = 0.0, solution_gap = 0.0;
      int evaluations = 0;
      auto certificate = [&](const VectorXd& theta, const VectorXd& u) {
        ++evaluations;
        if (!theta.allFinite() || !u.allFinite() || u.minCoeff() < low - 1e-12 || u.maxCoeff() > high + 1e-12) return(false);
        VectorXd bounded = u.cwiseMax(low).cwiseMin(high);
        VectorXd residual = response - design * theta;
        double primal = 0.5 * theta.squaredNorm();
        for (int j = 0; j < m; ++j) primal += residual[j] * (quantile - (residual[j] < 0 ? 1.0 : 0.0)) / m;
        const double dual_value = response.dot(bounded) - 0.5 * (design.transpose() * bounded).squaredNorm();
        const double gap = primal - dual_value;
        if (std::isfinite(gap) && gap >= -1e-9 && gap <= tolerance * std::max(1.0, std::abs(primal))) {
          solution_theta = theta;
          solution_dual = bounded;
          solution_loss = primal;
          solution_gap = gap;
          return(true);
        }
        return(false);
      };
      auto active_data = [&](const VectorXd& theta, VectorXd& u, std::vector<int>& active) {
        VectorXd residual = response - design * theta;
        active.clear();
        u = VectorXd::Zero(m);
        for (int j = 0; j < m; ++j) {
          const double roundoff = 128 * std::numeric_limits<double>::epsilon() *
            std::max(1.0, std::max(std::abs(response[j]), std::abs(design.row(j).dot(theta))));
          if (std::abs(residual[j]) <= roundoff) active.push_back(j);
          else u[j] = residual[j] < 0 ? low : high;
        }
      };
      auto recover = [&](const VectorXd& theta) {
        VectorXd u(m);
        std::vector<int> active;
        active_data(theta, u, active);
        if (certificate(theta, u)) return(true);
        VectorXd rhs = theta - design.transpose() * u;
        for (int a : active) {
          const double norm = design.row(a).squaredNorm();
          if (norm == 0) continue;
          VectorXd candidate = u;
          candidate[a] = design.row(a).dot(rhs) / norm;
          if (certificate(theta, candidate)) return(true);
        }
        for (std::size_t a = 0; a < active.size(); ++a) {
          for (std::size_t b = a + 1; b < active.size(); ++b) {
            MatrixXd columns(2, 2);
            columns.col(0) = design.row(active[a]).transpose();
            columns.col(1) = design.row(active[b]).transpose();
            Eigen::FullPivLU<MatrixXd> solver(columns);
            if (solver.rank() != 2) continue;
            std::vector<int> remaining;
            for (std::size_t c = 0; c < active.size(); ++c) if (c != a && c != b) remaining.push_back(active[c]);
            if (remaining.size() >= 24) Rcpp::stop("Too many coincident profile residuals for finite face enumeration.");
            const std::size_t patterns = std::size_t(1) << remaining.size();
            for (std::size_t pattern = 0; pattern < patterns; ++pattern) {
              VectorXd candidate = u;
              for (std::size_t c = 0; c < remaining.size(); ++c) candidate[remaining[c]] = (pattern >> c) & 1 ? high : low;
              VectorXd value = solver.solve(theta - design.transpose() * candidate);
              candidate[active[a]] = value[0];
              candidate[active[b]] = value[1];
              if (certificate(theta, candidate)) return(true);
            }
          }
        }
        if (!active.empty()) {
          int pivot = active[0];
          for (int a : active) if (design.row(a).squaredNorm() > design.row(pivot).squaredNorm()) pivot = a;
          const double norm = design.row(pivot).norm();
          if (norm > 0) {
            VectorXd axis = design.row(pivot).transpose() / norm;
            VectorXd candidate = u;
            double total_low = 0;
            for (int a : active) total_low += std::min(design.row(a).dot(axis) * low, design.row(a).dot(axis) * high);
            double remaining = axis.dot(rhs) - total_low;
            for (int a : active) {
              const double weight = design.row(a).dot(axis);
              if (weight == 0) continue;
              const double lower = std::min(weight * low, weight * high);
              const double upper = std::max(weight * low, weight * high);
              const double increment = std::max(0.0, std::min(upper - lower, remaining));
              candidate[a] = (lower + increment) / weight;
              remaining -= increment;
            }
            if (certificate(theta, candidate)) return(true);
          }
        }
        return(false);
      };
      VectorXd first(m);
      for (int j = 0; j < m; ++j) first[j] = response[j] < 0 ? low : high;
      bool solved = certificate(design.transpose() * first, first);
      for (int j = 0; j < m && !solved; ++j) {
        const double norm = design.row(j).squaredNorm();
        if (norm == 0) continue;
        VectorXd origin = response[j] * design.row(j).transpose() / norm;
        VectorXd tangent(2);
        tangent << -design(j, 1), design(j, 0);
        tangent /= std::sqrt(norm);
        VectorXd line_response = response - design * origin;
        VectorXd line_design = design * tangent;
        const double location = scalar_profile_root(line_response, line_design, 1.0, quantile);
        VectorXd theta = origin + location * tangent;
        solved = recover(theta);
        if (solved) break;
        VectorXd u(m);
        std::vector<int> active;
        active_data(theta, u, active);
        if (active.size() >= 24) Rcpp::stop("Too many coincident profile residuals for adjacent face enumeration.");
        const std::size_t patterns = std::size_t(1) << active.size();
        for (std::size_t pattern = 0; pattern < patterns && !solved; ++pattern) {
          VectorXd candidate = u;
          for (std::size_t a = 0; a < active.size(); ++a) candidate[active[a]] = (pattern >> a) & 1 ? high : low;
          solved = certificate(design.transpose() * candidate, candidate);
        }
      }
      if (!solved) Rcpp::stop("The two dimensional profile face search found no certified optimum.");
      gamma.row(i) = (transform * solution_theta).transpose();
      for (int j = 0; j < m; ++j) dual[row_order[begin + j]] = solution_dual[j];
      losses[i] = solution_loss;
      gaps[i] = solution_gap;
      iterations[i] = evaluations;
      continue;
    }
    MatrixXd transformed = factor.solve(Zi.transpose());
    VectorXd diagonal = (Zi.array() * transformed.transpose().array()).rowwise().sum();
    VectorXd u = VectorXd::Zero(m), nuisance = VectorXd::Zero(q);
    bool converged = false;
    for (int iteration = 1; iteration <= max_iterations; ++iteration) {
      for (int j = 0; j < m; ++j) {
        double candidate;
        if (diagonal[j] > 0) {
          candidate = u[j] + (response[j] - Zi.row(j).dot(nuisance)) / diagonal[j];
        } else {
          candidate = response[j] < 0 ? (quantile - 1.0) / m : quantile / m;
        }
        candidate = std::max((quantile - 1.0) / m, std::min(quantile / m, candidate));
        nuisance.noalias() += (candidate - u[j]) * transformed.col(j);
        u[j] = candidate;
      }
      VectorXd residual = response - Zi * nuisance;
      double check = 0.0;
      for (int j = 0; j < m; ++j) {
        check += residual[j] * (quantile - (residual[j] < 0 ? 1.0 : 0.0)) / m;
      }
      const double quadratic = nuisance.dot(Lambda * nuisance);
      const double primal = check + 0.5 * quadratic;
      const double gap = primal - response.dot(u) + 0.5 * quadratic;
      if (gap < -1e-9 || !std::isfinite(gap)) Rcpp::stop("Invalid cluster duality gap.");
      if (gap <= tolerance * std::max(1.0, std::abs(primal))) {
        gamma.row(i) = nuisance.transpose();
        for (int j = 0; j < m; ++j) dual[row_order[begin + j]] = u[j];
        losses[i] = primal;
        gaps[i] = gap;
        iterations[i] = iteration;
        converged = true;
        break;
      }
    }
    if (!converged) Rcpp::stop("Exact cluster profile reached its iteration limit.");
  }
  return Rcpp::List::create(Rcpp::_["gamma"] = gamma, Rcpp::_["dual"] = dual,
                            Rcpp::_["loss"] = losses, Rcpp::_["gap"] = gaps,
                            Rcpp::_["iterations"] = iterations);
}

// [[Rcpp::export]]
Rcpp::List joint_profile_pdal_cpp(
    const Eigen::Map<MatrixXd> X, const Eigen::Map<MatrixXd> Z,
    const Eigen::Map<VectorXd> outcome, Rcpp::IntegerVector cluster,
    const Eigen::Map<VectorXd> inverse_size, const Eigen::Map<MatrixXd> Lambda,
    const Eigen::Map<VectorXd> coordinate_penalty, double quantile,
    Eigen::VectorXd theta, Eigen::MatrixXd gamma, Eigen::VectorXd dual,
    int max_iterations = 30000, double gap_tolerance = 1e-8,
    double step_ratio = 1.0, double initial_step = 0.1) {
  const int N = X.rows(), p = X.cols(), q = Z.cols(), n = gamma.rows();
  if (N < 2 || p < 1 || q < 1 || n < 2 || Z.rows() != N || outcome.size() != N ||
      cluster.size() != N || inverse_size.size() != N || theta.size() != p ||
      gamma.cols() != q || dual.size() != N || coordinate_penalty.size() != p ||
      Lambda.rows() != q || Lambda.cols() != q || !(quantile > 0 && quantile < 1) ||
      !(gap_tolerance > 0) || !(step_ratio > 0) || !(initial_step > 0) || max_iterations < 1 ||
      !X.allFinite() || !Z.allFinite() || !outcome.allFinite() || !Lambda.allFinite() ||
      !theta.allFinite() || !gamma.allFinite() || !dual.allFinite() ||
      !coordinate_penalty.allFinite() || coordinate_penalty.minCoeff() < 0) {
    Rcpp::stop("Invalid joint profile input.");
  }
  Eigen::LLT<MatrixXd> lambda_factor(Lambda);
  if (lambda_factor.info() != Eigen::Success) Rcpp::stop("Lambda must be positive definite.");
  for (int r = 0; r < N; ++r) {
    if (cluster[r] < 0 || cluster[r] >= n || !(inverse_size[r] > 0)) {
      Rcpp::stop("Invalid cluster membership or weight.");
    }
    dual[r] = std::max((quantile - 1) * inverse_size[r],
                       std::min(quantile * inverse_size[r], dual[r]));
  }
  auto forward = [&](const VectorXd& value, const MatrixXd& nuisance) {
    VectorXd prediction = X * value;
    for (int r = 0; r < N; ++r) prediction[r] += Z.row(r).dot(nuisance.row(cluster[r]));
    return prediction;
  };
  auto adjoint_gamma = [&](const VectorXd& value) {
    MatrixXd result = MatrixXd::Zero(n, q);
    for (int r = 0; r < N; ++r) result.row(cluster[r]).noalias() += value[r] * Z.row(r);
    return result;
  };
  VectorXd prediction = forward(theta, gamma);
  VectorXd theta_gradient = X.transpose() * dual;
  MatrixXd gamma_gradient = adjoint_gamma(dual);
  VectorXd penalty = n * coordinate_penalty;
  double step = initial_step, extrapolation = 1.0, primal = 0.0, dual_value = 0.0;
  double gap = R_PosInf, beta_kkt = R_PosInf, gamma_kkt = R_PosInf;
  int total_backtracks = 0, accepted_iterations = 0;
  std::vector<double> trace_iteration, trace_gap, trace_objective, trace_step;
  for (int iteration = 1; iteration <= max_iterations; ++iteration) {
    if (iteration % 100 == 0) Rcpp::checkUserInterrupt();
    VectorXd theta_new(p);
    for (int j = 0; j < p; ++j) {
      theta_new[j] = threshold_value(theta[j] + step * theta_gradient[j], step * penalty[j]) /
        (1.0 + step);
    }
    Eigen::LLT<MatrixXd> proximal_factor(MatrixXd::Identity(q, q) + step * Lambda);
    MatrixXd gamma_new = proximal_factor.solve((gamma + step * gamma_gradient).transpose()).transpose();
    VectorXd prediction_new = forward(theta_new, gamma_new);
    double trial_step = std::min(10.0, step * std::sqrt(1.0 + extrapolation));
    VectorXd dual_new(N), theta_gradient_new(p);
    MatrixXd gamma_gradient_new(n, q);
    bool accepted = false;
    for (int backtrack = 0; backtrack < 60; ++backtrack) {
      const double ratio = trial_step / step;
      VectorXd prediction_extrapolated = prediction_new + ratio * (prediction_new - prediction);
      dual_new = dual + step_ratio * trial_step * (outcome - prediction_extrapolated);
      for (int r = 0; r < N; ++r) {
        dual_new[r] = std::max((quantile - 1) * inverse_size[r],
                               std::min(quantile * inverse_size[r], dual_new[r]));
      }
      theta_gradient_new.noalias() = X.transpose() * dual_new;
      gamma_gradient_new = adjoint_gamma(dual_new);
      const double lhs = step_ratio * trial_step * trial_step *
        ((theta_gradient_new - theta_gradient).squaredNorm() +
         (gamma_gradient_new - gamma_gradient).squaredNorm());
      const double rhs = 0.99 * 0.99 * (dual_new - dual).squaredNorm();
      if (lhs <= rhs + 1e-28) {
        accepted = true;
        total_backtracks += backtrack;
        break;
      }
      trial_step *= 0.5;
    }
    if (!accepted) Rcpp::stop("The primal dual linesearch did not terminate.");
    extrapolation = trial_step / step;
    step = trial_step;
    theta = theta_new;
    gamma = gamma_new;
    prediction = prediction_new;
    dual = dual_new;
    theta_gradient = theta_gradient_new;
    gamma_gradient = gamma_gradient_new;
    if (iteration % 10 == 0 || iteration == max_iterations) {
      primal = 0.5 * theta.squaredNorm() + (penalty.array() * theta.array().abs()).sum();
      for (int r = 0; r < N; ++r) {
        const double residual = outcome[r] - prediction[r];
        primal += inverse_size[r] * residual * (quantile - (residual < 0 ? 1.0 : 0.0));
      }
      primal += 0.5 * (gamma.array() * (gamma * Lambda).array()).sum();
      double conjugate = 0.0;
      beta_kkt = 0.0;
      for (int j = 0; j < p; ++j) {
        const double excess = threshold_value(theta_gradient[j], penalty[j]);
        conjugate += 0.5 * excess * excess;
        const double stationarity = theta[j] - theta_gradient[j];
        const double defect = theta[j] == 0.0 ? std::max(0.0, std::abs(stationarity) - penalty[j]) :
          std::abs(stationarity + (theta[j] > 0 ? penalty[j] : -penalty[j]));
        beta_kkt = std::max(beta_kkt, defect / n);
      }
      MatrixXd gamma_conjugate = lambda_factor.solve(gamma_gradient.transpose()).transpose();
      conjugate += 0.5 * (gamma_gradient.array() * gamma_conjugate.array()).sum();
      dual_value = outcome.dot(dual) - conjugate;
      gap = (primal - dual_value) / n;
      gamma_kkt = (gamma * Lambda - gamma_gradient).cwiseAbs().maxCoeff();
      if (!std::isfinite(gap) || gap < -1e-10 * std::max(1.0, std::abs(primal / n))) {
        Rcpp::stop("Invalid primal dual gap.");
      }
      trace_iteration.push_back(iteration);
      trace_gap.push_back(gap);
      trace_objective.push_back(primal / n);
      trace_step.push_back(step);
      if (gap <= gap_tolerance * std::max(1.0, std::abs(primal / n))) {
        accepted_iterations = iteration;
        break;
      }
    }
  }
  if (accepted_iterations == 0) {
    Rcpp::stop("Joint profile optimisation reached its iteration limit with gap %.6g.", gap);
  }
  return Rcpp::List::create(
    Rcpp::_["theta"] = theta, Rcpp::_["gamma"] = gamma, Rcpp::_["dual"] = dual,
    Rcpp::_["primal"] = primal / n, Rcpp::_["dual_value"] = dual_value / n,
    Rcpp::_["duality_gap"] = gap, Rcpp::_["beta_kkt"] = beta_kkt,
    Rcpp::_["gamma_kkt"] = gamma_kkt, Rcpp::_["iterations"] = accepted_iterations,
    Rcpp::_["backtracks"] = total_backtracks,
    Rcpp::_["history"] = Rcpp::DataFrame::create(
      Rcpp::_["iteration"] = trace_iteration, Rcpp::_["gap"] = trace_gap,
      Rcpp::_["objective"] = trace_objective, Rcpp::_["step"] = trace_step));
}
