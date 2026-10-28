#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Rcpp::List submission_clime_anchors_cpp(const Eigen::Map<Eigen::MatrixXd> H, double gamma) {
  const int p = H.rows();
  if (H.cols() != p || gamma <= 0 || gamma >= 1) Rcpp::stop("Invalid CLIME anchor input.");
  Eigen::MatrixXd anchors = Eigen::MatrixXd::Zero(p, p);
  Rcpp::LogicalVector valid(p, false);
  Rcpp::NumericVector ridge(p, NA_REAL), residual(p, NA_REAL);
  const double scale = std::max(H.trace(), 1e-8);
  const Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(p, p);
  int remaining = p;
  for (int exponent = 0; exponent <= 12 && remaining > 0; ++exponent) {
    const double lambda = scale * std::pow(10.0, -exponent);
    Eigen::MatrixXd system = H;
    system.diagonal().array() += lambda;
    Eigen::LDLT<Eigen::MatrixXd> decomposition(system);
    if (decomposition.info() != Eigen::Success || !decomposition.isPositive()) Rcpp::stop("The anchor ridge system is not positive definite.");
    Eigen::MatrixXd inverse = decomposition.solve(identity);
    Eigen::MatrixXd product = H * inverse;
    for (int k = 0; k < p; ++k) {
      if (valid[k] || product(k, k) <= 1e-12) continue;
      Eigen::VectorXd value = (1 - gamma) * inverse.col(k) / product(k, k);
      Eigen::VectorXd defect = H * value;
      defect[k] -= 1;
      const double error = defect.cwiseAbs().maxCoeff();
      if (std::isfinite(error) && error <= gamma + 1e-9) {
        anchors.col(k) = value;
        valid[k] = true;
        ridge[k] = lambda;
        residual[k] = error;
        --remaining;
      }
    }
  }
  return Rcpp::List::create(Rcpp::_["anchors"] = anchors, Rcpp::_["valid"] = valid,
    Rcpp::_["ridge"] = ridge, Rcpp::_["residual"] = residual);
}
