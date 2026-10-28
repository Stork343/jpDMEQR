#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;
using Eigen::VectorXd;

// [[Rcpp::export]]
Rcpp::List submission_sqr_fista_cpp(const Eigen::Map<MatrixXd> X,
    const Eigen::Map<VectorXd> y, double tau, double h, double lambda,
    const Eigen::Map<VectorXd> penalty, Eigen::VectorXd beta,
    int max_iterations, double beta_tolerance, double kkt_tolerance, int max_backtrack) {
  const int N = X.rows(), p = X.cols();
  if (y.size() != N || beta.size() != p || penalty.size() != p ||
      h <= 0 || tau <= 0 || tau >= 1 || lambda < 0) Rcpp::stop("Invalid Gaussian quantile programme.");
  auto smooth = [&](const VectorXd& b, VectorXd* gradient) {
    VectorXd residual = y - X * b;
    VectorXd score(N);
    double value = 0;
    for (int i = 0; i < N; ++i) {
      const double z = residual[i] / h;
      value += residual[i] * (tau - R::pnorm(-z, 0, 1, true, false)) + h * R::dnorm(z, 0, 1, false);
      score[i] = 1 - R::pnorm(z, 0, 1, true, false) - tau;
    }
    if (gradient != nullptr) *gradient = X.transpose() * score / N;
    return value / N;
  };
  auto regularizer = [&](const VectorXd& b) { return lambda * (penalty.array() * b.array().abs()).sum(); };
  auto residual_kkt = [&](const VectorXd& b, const VectorXd& g) {
    double answer = 0;
    for (int j = 0; j < p; ++j) {
      double value = std::abs(g[j]);
      if (penalty[j] > 0) {
        value = std::abs(b[j]) > 1e-10 ? std::abs(g[j] + lambda * penalty[j] * (b[j] > 0 ? 1 : -1)) :
          std::max(std::abs(g[j]) - lambda * penalty[j], 0.0);
      }
      answer = std::max(answer, value);
    }
    return answer;
  };
  VectorXd vector = VectorXd::Constant(p, 1.0 / std::sqrt(static_cast<double>(p)));
  double eigenvalue = 0, previous_eigenvalue = 0;
  for (int iteration = 0; iteration < 100; ++iteration) {
    VectorXd action = X.transpose() * (X * vector) / N;
    double norm = action.norm();
    if (!std::isfinite(norm) || norm == 0) Rcpp::stop("Unidentified Gaussian quantile design.");
    vector = action / norm;
    eigenvalue = (X * vector).squaredNorm() / N;
    if (std::abs(eigenvalue - previous_eigenvalue) <= 1e-8 * (1 + std::abs(eigenvalue))) break;
    previous_eigenvalue = eigenvalue;
  }
  double lipschitz = std::max(eigenvalue * R::dnorm(0, 0, 1, false) / h, 1e-6);
  VectorXd extrapolated = beta, gradient(p), next(p), final_gradient(p);
  double acceleration = 1, objective = smooth(beta, &final_gradient) + regularizer(beta);
  double kkt = residual_kkt(beta, final_gradient), change = R_PosInf;
  bool converged = false;
  int completed = 0;
  Rcpp::NumericMatrix history(max_iterations, 6);
  for (int iteration = 0; iteration < max_iterations; ++iteration) {
    double value = smooth(extrapolated, &gradient);
    bool accepted = false;
    int backtrack = 0;
    for (; backtrack <= max_backtrack; ++backtrack) {
      VectorXd proposal = extrapolated - gradient / lipschitz;
      for (int j = 0; j < p; ++j) {
        next[j] = (proposal[j] > 0 ? 1 : -1) * std::max(std::abs(proposal[j]) - lambda * penalty[j] / lipschitz, 0.0);
      }
      VectorXd delta = next - extrapolated;
      double candidate = smooth(next, nullptr);
      if (candidate <= value + gradient.dot(delta) + lipschitz * delta.squaredNorm() / 2 + 1e-12) {
        accepted = true;
        break;
      }
      lipschitz *= 2;
    }
    if (!accepted) Rcpp::stop("Gaussian quantile majorisation failed.");
    change = (next - beta).cwiseAbs().maxCoeff();
    objective = smooth(next, &final_gradient) + regularizer(next);
    kkt = residual_kkt(next, final_gradient);
    history(iteration, 0) = iteration + 1;
    history(iteration, 1) = objective;
    history(iteration, 2) = change;
    history(iteration, 3) = 1 / lipschitz;
    history(iteration, 4) = backtrack;
    history(iteration, 5) = kkt;
    completed = iteration + 1;
    if (change < beta_tolerance * std::max(1.0, next.cwiseAbs().maxCoeff()) && kkt <= kkt_tolerance) {
      beta = next;
      converged = true;
      break;
    }
    bool restart = (extrapolated - next).dot(next - beta) > 0;
    double next_acceleration = restart ? 1 : (1 + std::sqrt(1 + 4 * acceleration * acceleration)) / 2;
    if (restart) acceleration = 1;
    extrapolated = next + ((acceleration - 1) / next_acceleration) * (next - beta);
    beta = next;
    acceleration = next_acceleration;
  }
  Rcpp::NumericMatrix used(completed, 6);
  for (int i = 0; i < completed; ++i) for (int j = 0; j < 6; ++j) used(i, j) = history(i, j);
  return Rcpp::List::create(Rcpp::_["beta"] = beta, Rcpp::_["objective"] = objective,
    Rcpp::_["score"] = final_gradient, Rcpp::_["kkt_residual"] = kkt,
    Rcpp::_["beta_change"] = change, Rcpp::_["converged"] = converged,
    Rcpp::_["iterations"] = completed, Rcpp::_["history"] = used);
}
