#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp14)]]

// [[Rcpp::export]]
Rcpp::List submission_qgee_geometry_cpp(const Eigen::Map<Eigen::MatrixXd> X,
                                       const Eigen::Map<Eigen::VectorXd> score,
                                       const Eigen::Map<Eigen::VectorXd> density,
                                       const Eigen::Map<Eigen::VectorXd> f0,
                                       const Rcpp::IntegerVector sizes,
                                       const Rcpp::List inverses) {
  if (X.rows() != score.size() || X.rows() != density.size() || X.rows() != f0.size())
    Rcpp::stop("QGEE geometry dimensions differ.");
  Eigen::MatrixXd transformed(X.rows(), X.cols());
  Eigen::VectorXd transformed_score(X.rows());
  int start = 0;
  for (int i = 0; i < sizes.size(); ++i) {
    const int m = sizes[i];
    if (m == 0) continue;
    if (start + m > X.rows()) Rcpp::stop("QGEE cluster index exceeds the design.");
    const Eigen::MatrixXd inverse = Rcpp::as<Eigen::MatrixXd>(inverses[m - 1]);
    transformed.middleRows(start, m).noalias() =
      f0.segment(start, m).asDiagonal() * inverse *
      density.segment(start, m).asDiagonal() * X.middleRows(start, m);
    transformed_score.segment(start, m).noalias() =
      f0.segment(start, m).asDiagonal() * inverse * score.segment(start, m);
    start += m;
  }
  if (start != X.rows()) Rcpp::stop("QGEE cluster sizes do not cover the design.");
  Eigen::MatrixXd hessian = X.transpose() * transformed;
  Eigen::VectorXd gradient = X.transpose() * transformed_score;
  return Rcpp::List::create(Rcpp::Named("hessian") = hessian,
                            Rcpp::Named("score") = gradient);
}
