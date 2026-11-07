#include <RcppEigen.h>
#include <vector>
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp17)]]

using Eigen::MatrixXd;
using Eigen::VectorXd;

// [[Rcpp::export]]
Rcpp::List clinical_intercept_qp(
    const Eigen::Map<MatrixXd> X, const Eigen::Map<MatrixXd> Z,
    const Eigen::Map<VectorXd> outcome, Rcpp::IntegerVector cluster,
    const Eigen::Map<VectorXd> inverse_size, const Eigen::Map<MatrixXd> Lambda,
    const Eigen::Map<VectorXd> coordinate_penalty, double quantile,
    VectorXd theta, MatrixXd gamma, VectorXd dual,
    int max_iterations = 30000, double gap_tolerance = 1e-8) {
  const int N = X.rows(), d = X.cols(), n = gamma.rows();
  if (N < 1 || d < 1 || n < 1 || Z.rows() != N || Z.cols() != 1 ||
      Lambda.rows() != 1 || Lambda.cols() != 1 || !(Lambda(0, 0) > 0) ||
      gamma.cols() != 1 || theta.size() != d || outcome.size() != N ||
      cluster.size() != N || inverse_size.size() != N || dual.size() != N ||
      coordinate_penalty.size() != d || coordinate_penalty.cwiseAbs().maxCoeff() != 0 ||
      (Z.array() - 1).abs().maxCoeff() > 1e-12 || !X.allFinite() || !outcome.allFinite() ||
      !theta.allFinite() || !gamma.allFinite() || !dual.allFinite() ||
      !(quantile > 0 && quantile < 1) || !(gap_tolerance > 0) || max_iterations < 1)
    Rcpp::stop("Invalid intercept quadratic program.");
  std::vector<std::vector<int>> groups(n);
  std::vector<int> previous_pattern(N, 9);
  for (int j = 0; j < N; ++j) {
    if (cluster[j] < 0 || cluster[j] >= n) Rcpp::stop("Invalid patient index.");
    groups[cluster[j]].push_back(j);
  }
  MatrixXd means = MatrixXd::Zero(d, n);
  for (int i = 0; i < n; ++i) {
    if (groups[i].empty()) Rcpp::stop("An empty patient group is present.");
    for (int j : groups[i]) {
      if (std::abs(inverse_size[j] - 1.0 / groups[i].size()) > 1e-12)
        Rcpp::stop("Patient weights must sum to one.");
      means.col(i).noalias() += inverse_size[j] * X.row(j).transpose();
    }
  }
  const double lambda = Lambda(0, 0), rho = 1, relaxation = 1.6;
  MatrixXd weighted_design = X.array().colwise() * inverse_size.array();
  MatrixXd H = MatrixXd::Identity(d, d) + rho * X.transpose() * weighted_design;
  H.noalias() -= (rho * rho / (lambda + rho)) * means * means.transpose();
  Eigen::LLT<MatrixXd> factor(H);
  if (factor.info() != Eigen::Success) Rcpp::stop("The quadratic system is not positive definite.");
  VectorXd effect = gamma.col(0);
  VectorXd lower = (quantile-1) * inverse_size, upper = quantile * inverse_size;
  dual = dual.cwiseMax(lower).cwiseMin(upper);
  VectorXd u = -dual.array() / (rho * inverse_size.array());
  VectorXd prediction = X * theta;
  for (int j = 0; j < N; ++j) prediction[j] += effect[cluster[j]];
  VectorXd residual = outcome - prediction;
  double primal = R_PosInf, dual_value = R_NegInf, gap = R_PosInf;
  int face_checks = 0;
  auto evaluate = [&](const VectorXd& coefficient, const VectorXd& nuisance,
                      const VectorXd& score, double& value, double& bound) {
    VectorXd r = outcome - X * coefficient;
    VectorXd group_score = VectorXd::Zero(n);
    double loss = 0;
    for (int j = 0; j < N; ++j) {
      r[j] -= nuisance[cluster[j]];
      loss += inverse_size[j] * r[j] * (quantile - (r[j] < 0));
      group_score[cluster[j]] += score[j];
    }
    value = (loss + coefficient.squaredNorm()/2 + lambda*nuisance.squaredNorm()/2)/n;
    bound = (outcome.dot(score) - (X.transpose()*score).squaredNorm()/2 -
             group_score.squaredNorm()/(2*lambda))/n;
    return value-bound;
  };
  auto finish = [&](int iterations, bool polished) {
    VectorXd group_score = VectorXd::Zero(n);
    for (int j = 0; j < N; ++j) group_score[cluster[j]] += dual[j];
    gamma.col(0) = effect;
    return Rcpp::List::create(Rcpp::_["theta"]=theta,Rcpp::_["gamma"]=gamma,Rcpp::_["dual"]=dual,
      Rcpp::_["primal"]=primal,Rcpp::_["dual_value"]=dual_value,Rcpp::_["duality_gap"]=gap,
      Rcpp::_["beta_kkt"]=(theta-X.transpose()*dual).cwiseAbs().maxCoeff()/n,
      Rcpp::_["gamma_kkt"]=(lambda*effect-group_score).cwiseAbs().maxCoeff(),
      Rcpp::_["iterations"]=iterations,Rcpp::_["face_checks"]=face_checks,Rcpp::_["polished"]=polished,
      Rcpp::_["solver"]="bounded_intercept_qp",Rcpp::_["solver_status"]="certified");
  };
  for (int iteration = 1; iteration <= max_iterations; ++iteration) {
    Rcpp::checkUserInterrupt();
    VectorXd v = outcome-residual-u, mean_v = VectorXd::Zero(n);
    for (int j = 0; j < N; ++j) mean_v[cluster[j]] += inverse_size[j]*v[j];
    VectorXd right = rho * weighted_design.transpose()*v -
                     (rho*rho/(lambda+rho))*means*mean_v;
    theta = factor.solve(right);
    effect = (rho/(lambda+rho))*(mean_v-means.transpose()*theta);
    prediction.noalias() = X*theta;
    for (int j = 0; j < N; ++j) prediction[j] += effect[cluster[j]];
    VectorXd relaxed = relaxation*prediction+(1-relaxation)*(outcome-residual);
    VectorXd argument = outcome-relaxed-u;
    for (int j = 0; j < N; ++j) {
      if (argument[j] > quantile/rho) residual[j] = argument[j]-quantile/rho;
      else if (argument[j] < (quantile-1)/rho) residual[j] = argument[j]-(quantile-1)/rho;
      else residual[j] = 0;
    }
    u += relaxed+residual-outcome;
    dual = (-rho*inverse_size.array()*u.array()).matrix().cwiseMax(lower).cwiseMin(upper);
    gap = evaluate(theta,effect,dual,primal,dual_value);
    if (!std::isfinite(gap) || gap < -1e-10*std::max(1.0,std::abs(primal)))
      Rcpp::stop("The quadratic certificate is invalid.");
    if (gap <= gap_tolerance*std::max(1.0,std::abs(primal))) return finish(iteration,false);
    if (iteration != 1 && iteration%5 != 0) continue;
    std::vector<int> pattern(N), anchors(n,-1), active, constrained;
    bool changed = false;
    for (int j = 0; j < N; ++j) {
      pattern[j] = residual[j]>0 ? 1 : (residual[j]<0 ? -1 : 0);
      changed = changed || pattern[j]!=previous_pattern[j];
    }
    if (!changed) continue;
    previous_pattern = pattern;
    right.setZero();
    for (int i = 0; i < n; ++i) {
      for (int j : groups[i]) {
        if (pattern[j]==0 && anchors[i]<0) {
          anchors[i]=j;
          active.push_back(j);
        } else if (pattern[j]==0) constrained.push_back(j);
      }
      const int anchor = anchors[i];
      if (anchor>=0) right.noalias() += lambda*outcome[anchor]*X.row(anchor).transpose();
      for (int j : groups[i]) {
        if (pattern[j]==0) continue;
        double score = pattern[j]>0 ? upper[j] : lower[j];
        right.noalias() += score*X.row(j).transpose();
        if (anchor>=0) right.noalias() -= score*X.row(anchor).transpose();
      }
    }
    if (static_cast<int>(constrained.size())>d) continue;
    ++face_checks;
    MatrixXd active_design(active.size(),d);
    for (unsigned int j = 0; j < active.size(); ++j) active_design.row(j)=X.row(active[j]);
    MatrixXd face_H = MatrixXd::Identity(d,d)+lambda*active_design.transpose()*active_design;
    Eigen::LLT<MatrixXd> face_factor(face_H);
    if (face_factor.info()!=Eigen::Success) Rcpp::stop("The active quadratic system is invalid.");
    VectorXd candidate=face_factor.solve(right), multipliers(constrained.size());
    if (!constrained.empty()) {
      MatrixXd constraints(constrained.size(),d);
      VectorXd values(constrained.size());
      for (unsigned int k=0; k<constrained.size(); ++k) {
        const int j=constrained[k], anchor=anchors[cluster[j]];
        constraints.row(k)=X.row(anchor)-X.row(j);
        values[k]=outcome[anchor]-outcome[j];
      }
      MatrixXd inverse_constraints=face_factor.solve(constraints.transpose());
      MatrixXd schur=constraints*inverse_constraints;
      Eigen::ColPivHouseholderQR<MatrixXd> decomposition(schur);
      decomposition.setThreshold(1e-12);
      if (decomposition.rank()!=schur.rows()) continue;
      multipliers=decomposition.solve(constraints*candidate-values);
      candidate.noalias()-=inverse_constraints*multipliers;
    }
    VectorXd candidate_effect(n), candidate_dual=VectorXd::Zero(N);
    for (int j=0; j<N; ++j)
      if (pattern[j]!=0) candidate_dual[j]=pattern[j]>0 ? upper[j] : lower[j];
    for (unsigned int k=0; k<constrained.size(); ++k) candidate_dual[constrained[k]]=multipliers[k];
    for (int i=0; i<n; ++i) {
      double total=0;
      for (int j : groups[i]) total+=candidate_dual[j];
      const int anchor=anchors[i];
      if (anchor>=0) {
        candidate_effect[i]=outcome[anchor]-X.row(anchor).dot(candidate);
        candidate_dual[anchor]=lambda*candidate_effect[i]-total;
      } else candidate_effect[i]=total/lambda;
    }
    candidate_dual=candidate_dual.cwiseMax(lower).cwiseMin(upper);
    double candidate_primal, candidate_bound;
    const double candidate_gap=evaluate(candidate,candidate_effect,candidate_dual,candidate_primal,candidate_bound);
    if (std::isfinite(candidate_gap) && candidate_gap>=-1e-10*std::max(1.0,std::abs(candidate_primal)) &&
        candidate_gap<=gap_tolerance*std::max(1.0,std::abs(candidate_primal))) {
      theta.swap(candidate);
      effect.swap(candidate_effect);
      dual.swap(candidate_dual);
      primal=candidate_primal;
      dual_value=candidate_bound;
      gap=candidate_gap;
      return finish(iteration,true);
    }
  }
  Rcpp::stop("The intercept solver exhausted its certified iteration budget.");
}
