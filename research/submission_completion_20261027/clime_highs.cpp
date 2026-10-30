#include <RcppEigen.h>
#include <Highs.h>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <limits>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp17)]]

namespace {
void require_highs(HighsStatus status) {
  if (status == HighsStatus::kError) Rcpp::stop("HiGHS rejected the CLIME model or numerical control.");
}

std::vector<int> violated_coordinates(const Eigen::VectorXd& value,
                                     const std::vector<bool>& selected,
                                     double threshold) {
  std::vector<int> result;
  for (int i = 0; i < value.size(); ++i)
    if (!selected[i] && std::abs(value[i]) > threshold) result.push_back(i);
  std::stable_sort(result.begin(), result.end(), [&](int left, int right) {
    return std::abs(value[left]) > std::abs(value[right]);
  });
  if (result.size() > 16) result.resize(16);
  return result;
}

void append_coordinates(std::vector<int>& indices, std::vector<bool>& present,
                        const std::vector<int>& additions) {
  for (int value : additions) if (!present[value]) {
    indices.push_back(value);
    present[value] = true;
  }
}

double box_loss_bound(const Eigen::VectorXd& residual, const Eigen::MatrixXd& columns,
                      const Eigen::VectorXd& radii) {
  if (!columns.cols()) return residual.squaredNorm();
  auto evaluate = [&](const Eigen::VectorXd& direction) {
    const double norm = direction.squaredNorm();
    if (norm == 0) return 0.0;
    if (!std::isfinite(norm)) Rcpp::stop("The validation certificate direction is nonfinite.");
    const Eigen::VectorXd products = columns.transpose() * direction;
    const double support = radii.dot(products.cwiseAbs());
    const double linear = residual.dot(direction);
    const double magnitude = residual.cwiseAbs().dot(direction.cwiseAbs()) +
      radii.dot(columns.cwiseAbs().transpose() * direction.cwiseAbs());
    if (!std::isfinite(support) || !std::isfinite(linear) || !std::isfinite(magnitude))
      Rcpp::stop("The validation certificate arithmetic is nonfinite.");
    const double rounding = 64 * std::numeric_limits<double>::epsilon() *
      (residual.size() + radii.size() + 1) * magnitude;
    const double value = std::max(0.0, linear - support - rounding);
    const double denominator = norm * (1 + 64 * std::numeric_limits<double>::epsilon() * (residual.size() + 1));
    const double bound = value * value / denominator;
    if (!std::isfinite(bound)) Rcpp::stop("The validation loss bound is nonfinite.");
    return bound;
  };
  const Eigen::VectorXd coefficients = columns.colPivHouseholderQr().solve(-residual);
  const Eigen::VectorXd projected = residual + columns * coefficients;
  if (!projected.allFinite()) Rcpp::stop("The validation projection is nonfinite.");
  return std::max(evaluate(residual), evaluate(projected));
}
}

// [[Rcpp::export]]
Rcpp::List submission_clime_highs_cpp(const Eigen::Map<Eigen::MatrixXd> H,
                                     double gamma, bool presolve,
                                     Rcpp::Nullable<Rcpp::NumericMatrix> validation,
                                     double loss_limit) {
  const int p = H.rows();
  if (p < 1 || H.cols() != p || !H.allFinite() || !std::isfinite(gamma) || gamma <= 0 || gamma >= 1)
    Rcpp::stop("Invalid CLIME Hessian or tolerance.");
  const double scale = std::max(1.0, H.cwiseAbs().maxCoeff());
  if ((H - H.transpose()).cwiseAbs().maxCoeff() > 1e-10 * scale)
    Rcpp::stop("The CLIME Hessian must be symmetric.");

  Highs solver;
  require_highs(solver.setOptionValue("output_flag", false));
  require_highs(solver.setOptionValue("threads", 1));
  require_highs(solver.setOptionValue("solver", "simplex"));
  require_highs(solver.setOptionValue("presolve", presolve ? "on" : "off"));
  require_highs(solver.setOptionValue("simplex_strategy", 1));
  require_highs(solver.setOptionValue("simplex_scale_strategy", 2));
  require_highs(solver.setOptionValue("simplex_iteration_limit", 20 * p));
  require_highs(solver.setOptionValue("primal_feasibility_tolerance", 1e-8));
  require_highs(solver.setOptionValue("dual_feasibility_tolerance", 1e-8));

  Eigen::MatrixXd raw = Eigen::MatrixXd::Constant(p, p, NA_REAL);
  Eigen::MatrixXd duals = Eigen::MatrixXd::Constant(p, p, NA_REAL);
  Rcpp::NumericVector residuals(p, NA_REAL), dual_residuals(p, NA_REAL), gaps(p, NA_REAL);
  Rcpp::NumericVector objectives(p, NA_REAL), iterations(p, 0.0);
  Rcpp::IntegerVector programmes(p, 0), row_counts(p, 0), column_counts(p, 0);
  Rcpp::CharacterVector statuses(p, "not_required_after_global_infeasibility");
  int failed_coordinate = 0;
  double ray_residual = NA_REAL, ray_margin = NA_REAL;
  std::vector<bool> solved(p, false);
  std::vector<std::vector<int>> supports(p);
  std::vector<int> solved_order;

  auto solve_coordinate = [&](int k) {
    if (solved[k] || failed_coordinate) return;
    Rcpp::checkUserInterrupt();
    Eigen::VectorXd omega = Eigen::VectorXd::Zero(p);
    Eigen::VectorXd dual = Eigen::VectorXd::Zero(p);
    Eigen::VectorXd defect(p), pricing(p);
    Eigen::Index largest;
    const double largest_value = H.row(k).cwiseAbs().maxCoeff(&largest);
    bool certified = false;
    if (largest_value > 1e-14) {
      omega[largest] = (1 - gamma) / H(k, largest);
      dual[k] = 1 / largest_value;
      defect = H.col(largest) * omega[largest];
      defect[k] -= 1;
      pricing = H.col(k) * dual[k];
      certified = defect.cwiseAbs().maxCoeff() <= gamma + 1e-7;
    }

    std::vector<int> rows(1, k), columns(1, k);
    std::vector<bool> row_present(p, false), column_present(p, false);
    row_present[k] = column_present[k] = true;

    for (int step = 0; !certified && step < 2 * p; ++step) {
      Rcpp::checkUserInterrupt();
      const int r = rows.size(), d = columns.size();
      HighsLp model;
      model.num_col_ = 2 * d;
      model.num_row_ = r;
      model.sense_ = ObjSense::kMinimize;
      model.col_cost_.assign(2 * d, 1.0);
      model.col_lower_.assign(2 * d, 0.0);
      model.col_upper_.assign(2 * d, kHighsInf);
      model.row_lower_.resize(r);
      model.row_upper_.resize(r);
      for (int i = 0; i < r; ++i) {
        const double target = rows[i] == k ? 1.0 : 0.0;
        model.row_lower_[i] = target - gamma;
        model.row_upper_[i] = target + gamma;
      }
      model.a_matrix_.format_ = MatrixFormat::kColwise;
      model.a_matrix_.start_.clear();
      model.a_matrix_.start_.reserve(2 * d + 1);
      model.a_matrix_.index_.reserve(2 * d * r);
      model.a_matrix_.value_.reserve(2 * d * r);
      for (int j = 0; j < 2 * d; ++j) {
        model.a_matrix_.start_.push_back(model.a_matrix_.value_.size());
        const int original = columns[j % d];
        const double sign = j < d ? 1.0 : -1.0;
        for (int i = 0; i < r; ++i) {
          const double value = sign * H(rows[i], original);
          if (value != 0.0) {
            model.a_matrix_.index_.push_back(i);
            model.a_matrix_.value_.push_back(value);
          }
        }
      }
      model.a_matrix_.start_.push_back(model.a_matrix_.value_.size());
      model.setMatrixDimensions();
      require_highs(solver.passModel(std::move(model)));
      require_highs(solver.run());
      ++programmes[k];
      iterations[k] += solver.getInfo().simplex_iteration_count;
      row_counts[k] = r;
      column_counts[k] = d;
      const HighsModelStatus status = solver.getModelStatus();
      if (status == HighsModelStatus::kInfeasible) {
        bool has_ray = false;
        std::vector<double> reduced_ray(r);
        require_highs(solver.getDualRay(has_ray, reduced_ray.data()));
        if (!has_ray) Rcpp::stop("HiGHS did not provide a CLIME infeasibility certificate.");
        Eigen::VectorXd ray = Eigen::VectorXd::Zero(p);
        for (int i = 0; i < r; ++i) ray[rows[i]] = reduced_ray[i];
        const double norm = ray.cwiseAbs().sum();
        if (!ray.allFinite() || norm <= 0) Rcpp::stop("The CLIME dual ray is invalid.");
        ray /= norm;
        if (ray[k] < 0) ray = -ray;
        ray_residual = (H * ray).cwiseAbs().maxCoeff();
        ray_margin = ray[k] - gamma;
        if (ray_residual > 1e-8 * scale || ray_margin <= 1e-7)
          Rcpp::stop("The CLIME infeasibility certificate fails the complete matrix conditions.");
        duals.row(k) = ray.transpose();
        statuses[k] = "infeasible";
        failed_coordinate = k + 1;
        break;
      }
      if (status != HighsModelStatus::kOptimal)
        Rcpp::stop("The compiled CLIME subproblem did not attain an accepted optimum.");
      const HighsSolution& solution = solver.getSolution();
      if (solution.col_value.size() != static_cast<size_t>(2 * d) ||
          solution.row_dual.size() != static_cast<size_t>(r))
        Rcpp::stop("The CLIME solution dimensions differ.");
      omega.setZero();
      dual.setZero();
      defect.setZero();
      defect[k] = -1;
      pricing.setZero();
      for (int j = 0; j < d; ++j) {
        const double value = solution.col_value[j] - solution.col_value[d + j];
        omega[columns[j]] = value;
        if (value != 0.0) defect.noalias() += value * H.col(columns[j]);
      }
      for (int i = 0; i < r; ++i) {
        const double value = solution.row_dual[i];
        dual[rows[i]] = value;
        if (value != 0.0) pricing.noalias() += value * H.col(rows[i]);
      }
      if (!omega.allFinite() || !dual.allFinite() || !defect.allFinite() || !pricing.allFinite())
        Rcpp::stop("The compiled CLIME solution is nonfinite.");
      const double objective = omega.cwiseAbs().sum();
      const double gap = objective - dual[k] + gamma * dual.cwiseAbs().sum();
      certified = defect.cwiseAbs().maxCoeff() <= gamma + 1e-7 &&
        pricing.cwiseAbs().maxCoeff() <= 1 + 1e-7 &&
        std::abs(gap) <= 1e-7 * std::max(1.0, objective);
      if (certified) break;
      const std::vector<int> more_rows = violated_coordinates(defect, row_present, gamma + 1e-7);
      const std::vector<int> more_columns = violated_coordinates(pricing, column_present, 1 + 1e-7);
      if (more_rows.empty() && more_columns.empty())
        Rcpp::stop("CLIME certificate at coordinate %d and gamma %.10g: primal %.10g, dual %.10g, gap %.10g, objective %.10g, rows %d, columns %d.",
          k + 1, gamma, defect.cwiseAbs().maxCoeff(), pricing.cwiseAbs().maxCoeff(), gap, objective, r, d);
      append_coordinates(rows, row_present, more_rows);
      append_coordinates(columns, column_present, more_rows);
      append_coordinates(columns, column_present, more_columns);
    }
    if (failed_coordinate) return;
    if (!certified) Rcpp::stop("The compiled CLIME working sets exhausted their refinement budget.");
    raw.row(k) = omega.transpose();
    duals.row(k) = dual.transpose();
    objectives[k] = omega.cwiseAbs().sum();
    residuals[k] = defect.cwiseAbs().maxCoeff();
    dual_residuals[k] = pricing.cwiseAbs().maxCoeff();
    gaps[k] = objectives[k] - dual[k] + gamma * dual.cwiseAbs().sum();
    if (dual_residuals[k] > 1 + 1e-7 || std::abs(gaps[k]) > 1e-7 * std::max(1.0, objectives[k]))
      Rcpp::stop("The CLIME analytic direction fails its full dual certificate.");
    statuses[k] = "optimal";
    solved[k] = true;
    solved_order.push_back(k);
    for (int j = 0; j < p; ++j) if (omega[j] != 0.0) supports[k].push_back(j);
  };

  Rcpp::NumericVector column_losses(p, NA_REAL);
  Rcpp::NumericVector column_bounds(p, 0.0);
  int evaluated_columns = 0;
  long double lower_bound = 0;
  bool pruned = false;
  if (validation.isNotNull()) {
    Rcpp::NumericMatrix validation_input(validation.get());
    if (validation_input.nrow() != p || validation_input.ncol() != p ||
        std::isnan(loss_limit) || loss_limit < 0)
      Rcpp::stop("Invalid CLIME validation design or loss bound.");
    const Eigen::Map<Eigen::MatrixXd> V(validation_input.begin(), p, p);
    if (!V.allFinite()) Rcpp::stop("The CLIME validation Hessian is nonfinite.");
    const Eigen::VectorXd validation_norms = V.colwise().norm().transpose();
    std::vector<int> last_unknown(p, -1);
    auto update_bounds = [&]() {
      for (int k : solved_order) {
        if (!Rcpp::NumericVector::is_na(column_losses[k])) continue;
        std::vector<int> unknown;
        for (int j : supports[k]) if (!solved[j]) unknown.push_back(j);
        if (last_unknown[k] >= 0 && !unknown.empty() &&
            4 * static_cast<int>(unknown.size()) > 3 * last_unknown[k]) continue;
        Eigen::VectorXd residual = Eigen::VectorXd::Zero(p);
        residual[k] = -1;
        for (int j : supports[k]) if (solved[j]) {
          const int first = std::min(j, k), second = std::max(j, k);
          const double value = std::abs(raw(first, second)) <= std::abs(raw(second, first)) ?
            raw(first, second) : raw(second, first);
          if (value != 0.0) residual.noalias() += value * V.col(j);
        }
        double bound;
        if (unknown.empty()) {
          bound = residual.squaredNorm() / p;
          column_losses[k] = bound;
          ++evaluated_columns;
          if (column_bounds[k] > bound + 1e-10 * std::max(1.0, bound))
            Rcpp::stop("A validation certificate exceeds its completed column loss.");
        } else {
          Eigen::MatrixXd remaining(p, unknown.size());
          Eigen::VectorXd radii(unknown.size());
          for (int j = 0; j < static_cast<int>(unknown.size()); ++j) {
            remaining.col(j) = V.col(unknown[j]);
            radii[j] = std::abs(raw(k, unknown[j]));
          }
          bound = std::max(static_cast<double>(column_bounds[k]), box_loss_bound(residual, remaining, radii) / p);
        }
        if (!std::isfinite(bound) || bound < 0) Rcpp::stop("A validation column bound is invalid.");
        lower_bound += bound - column_bounds[k];
        column_bounds[k] = bound;
        last_unknown[k] = unknown.size();
        if (std::isfinite(loss_limit) && lower_bound > loss_limit + 1e-8 * std::max(1.0, loss_limit)) {
          pruned = true;
          return;
        }
      }
    };
    std::vector<int> order(p);
    std::iota(order.begin(), order.end(), 0);
    Eigen::VectorXd priority(p);
    for (int k = 0; k < p; ++k) {
      const double coefficient = (1 - gamma) / std::max(std::abs(H(k, k)), 1e-14);
      priority[k] = coefficient * coefficient * validation_norms[k] * validation_norms[k] - 2 * coefficient * V(k, k) + 1;
    }
    std::stable_sort(order.begin(), order.end(), [&](int left, int right) {
      return priority[left] > priority[right];
    });
    if (!std::isfinite(loss_limit)) {
      for (int k = 0; k < p; ++k) {
        solve_coordinate(k);
        if (failed_coordinate) break;
      }
      if (!failed_coordinate) update_bounds();
    } else for (int k : order) {
      if (!Rcpp::NumericVector::is_na(column_losses[k])) continue;
      solve_coordinate(k);
      if (failed_coordinate) break;
      update_bounds();
      if (pruned) break;
      std::vector<int> neighbours = supports[k];
      std::stable_sort(neighbours.begin(), neighbours.end(), [&](int left, int right) {
        return std::abs(raw(k, left)) * validation_norms[left] > std::abs(raw(k, right)) * validation_norms[right];
      });
      for (int j : neighbours) if (!solved[j]) {
        solve_coordinate(j);
        if (failed_coordinate) break;
        update_bounds();
        if (pruned) break;
      }
      if (failed_coordinate || pruned) break;
    }
  } else {
    for (int k = 0; k < p; ++k) {
      solve_coordinate(k);
      if (failed_coordinate) break;
    }
  }
  const int solved_rows = std::count(solved.begin(), solved.end(), true);
  if (pruned) for (int k = 0; k < p; ++k)
    if (!solved[k]) statuses[k] = "not_required_after_loss_bound";
  return Rcpp::List::create(
    Rcpp::_["omega_raw"] = pruned ? R_NilValue : Rcpp::wrap(raw),
    Rcpp::_["dual"] = pruned ? R_NilValue : Rcpp::wrap(duals),
    Rcpp::_["feasible"] = failed_coordinate == 0 && solved_rows == p,
    Rcpp::_["failed_coordinate"] = failed_coordinate,
    Rcpp::_["status"] = statuses, Rcpp::_["residual"] = residuals,
    Rcpp::_["dual_residual"] = dual_residuals, Rcpp::_["dual_gap"] = gaps,
    Rcpp::_["objective"] = objectives, Rcpp::_["programmes"] = programmes,
    Rcpp::_["simplex_iterations"] = iterations, Rcpp::_["rows"] = row_counts,
    Rcpp::_["columns"] = column_counts, Rcpp::_["ray_residual"] = ray_residual,
    Rcpp::_["ray_margin"] = ray_margin, Rcpp::_["highs_version"] = solver.version(),
    Rcpp::_["pruned"] = pruned, Rcpp::_["solved_rows"] = solved_rows,
    Rcpp::_["evaluated_columns"] = evaluated_columns, Rcpp::_["column_losses"] = column_losses,
    Rcpp::_["column_lower_bounds"] = column_bounds,
    Rcpp::_["loss_lower_bound"] = validation.isNotNull() ? static_cast<double>(lower_bound) : NA_REAL);
}

// [[Rcpp::export]]
double submission_sqr_box_bound_cpp(const Eigen::Map<Eigen::VectorXd> residual,
                                   const Eigen::Map<Eigen::MatrixXd> columns,
                                   const Eigen::Map<Eigen::VectorXd> radii) {
  if (columns.rows() != residual.size() || columns.cols() != radii.size() ||
      !residual.allFinite() || !columns.allFinite() || !radii.allFinite() || (radii.array() < 0).any())
    Rcpp::stop("Invalid validation box certificate inputs.");
  return box_loss_bound(residual, columns, radii);
}

// [[Rcpp::export]]
double submission_sqr_sparse_loss_cpp(const Eigen::Map<Eigen::MatrixXd> H,
                                      const Eigen::Map<Eigen::MatrixXd> W) {
  const int p = H.rows();
  if (H.cols() != p || W.rows() != p || W.cols() != p || !H.allFinite() || !W.allFinite())
    Rcpp::stop("Invalid SQR validation matrices.");
  Eigen::SparseMatrix<double> sparse = W.sparseView(0.0, 0.0);
  Eigen::MatrixXd defect = H * sparse;
  defect.diagonal().array() -= 1;
  return defect.squaredNorm() / p;
}
