/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_SPL_SPLINE_H_
#define SRC_SPL_SPLINE_H_

#include <algorithm>
#include <array>
#include <functional>
#include <numeric>
#include <sstream>
#include <utility>
#include <vector>

#include "control_point.h"
#include "knot_vector.h"
#include "multi_index_handler.h"
#include "parameter_space.h"
#include "physical_space.h"

namespace spl {
template<int DIM>
class Spline {
 public:
  virtual ~Spline() = default;
  Spline() = default;
  Spline(KnotVectors<DIM> knot_vector, std::array<Degree, DIM> degree) {
    parameter_space_ = std::make_shared<ParameterSpace<DIM>>(ParameterSpace<DIM>(knot_vector, degree));
  }
  explicit Spline(std::shared_ptr<ParameterSpace<DIM>> parameter_space) {
    parameter_space_ = parameter_space;
  }
  Spline(const Spline<DIM> &spline) {
    ParameterSpace<DIM> parameter_space(*spline.parameter_space_);
    parameter_space_ = std::make_shared<ParameterSpace<DIM>>(parameter_space);
  }

  virtual std::vector<double> Evaluate(std::array<ParamCoord, DIM> param_coord,
                                       const std::vector<int> &dimensions) const {
    this->ThrowIfParametricCoordinateOutsideKnotVectorRange(param_coord);

    auto first_non_zero = GetArrayOfFirstNonZeroBasisFunctions(param_coord);
    util::MultiIndexHandler<DIM> basisFunctionHandler(this->GetNumberOfBasisFunctionsToEvaluate());
    std::vector<double> evaluated_point(dimensions.size(), 0);

    for (int i = 0; i < basisFunctionHandler.Get1DLength(); ++i, basisFunctionHandler++) {
      auto indices = basisFunctionHandler.GetIndices();
      std::transform(indices.begin(), indices.end(), first_non_zero.begin(), indices.begin(), std::plus<>());
      for (size_t j = 0; j < dimensions.size(); ++j) {
        evaluated_point[j] += GetEvaluatedControlPoint(param_coord, indices, dimensions[j]);
      }
    }
    return evaluated_point;
  }

  virtual std::vector<double> EvaluateDerivative(std::array<ParamCoord, DIM> param_coord,
                                                 const std::vector<int> &dimensions,
                                                 std::array<int, DIM> derivative) const {
    this->ThrowIfParametricCoordinateOutsideKnotVectorRange(param_coord);

    auto first_non_zero = this->GetArrayOfFirstNonZeroBasisFunctions(param_coord);
    util::MultiIndexHandler<DIM> basisFunctionHandler(this->GetNumberOfBasisFunctionsToEvaluate());
    std::vector<double> evaluated_point(dimensions.size(), 0);

    for (int i = 0; i < basisFunctionHandler.Get1DLength(); ++i, basisFunctionHandler++) {
      auto indices = basisFunctionHandler.GetIndices();
      std::transform(indices.begin(), indices.end(), first_non_zero.begin(), indices.begin(), std::plus<>());
      for (size_t j = 0; j < dimensions.size(); ++j) {
        evaluated_point[j] += GetEvaluatedDerivativeControlPoint(param_coord, derivative, indices, dimensions[j]);
      }
    }
    return evaluated_point;
  }

  [[nodiscard]] std::vector<double> EvaluateAllNonZeroBasisFunctions(int direction, ParamCoord param_coord) const {
    return parameter_space_->EvaluateAllNonZeroBasisFunctions(direction, param_coord);
  }

  [[nodiscard]] std::vector<double> EvaluateAllNonZeroBasisFunctionDerivatives(int direction,
                                                                               ParamCoord param_coord,
                                                                               int derivative) const {
    return parameter_space_->EvaluateAllNonZeroBasisFunctionDerivatives(direction, param_coord, derivative);
  }

  bool AreGeometricallyEqual(const spl::Spline<DIM> &rhs,
                             double tolerance = util::NumericSettings<double>::kEpsilon()) const {
    double number = ceil(pow(100, 1.0 / DIM));
    std::array<int, DIM> points = std::array<int, DIM>();
    std::vector<int> dimensions;
    for (int i = 0; i < DIM; ++i) {
      points[i] = number + 1;
    }
    for (int i = 0; i < GetPointDim(); ++i) {
      dimensions.push_back(i);
    }
    util::MultiIndexHandler<DIM> point_handler(points);
    for (int i = 0; i < point_handler.Get1DLength(); ++i, ++point_handler) {
      std::array<ParamCoord, DIM> param_coord;
      for (int dim = 0; dim < DIM; ++dim) {
        double span = GetKnotVector(dim)->GetLastKnot().get() - GetKnotVector(dim)->GetKnot(0).get();
        param_coord[dim] = ParamCoord(span / number * point_handler[dim]);
      }
      std::vector<double> evaluate_this = Evaluate(param_coord, dimensions);
      std::vector<double> evaluate_rhs = rhs.Evaluate(param_coord, dimensions);
      if (util::VectorUtils<double>::ComputeDistance(evaluate_this, evaluate_rhs) > tolerance) {
        return false;
      }
    }
    return true;
  }

  [[nodiscard]] Degree GetDegree(int i) const {
    return parameter_space_->GetDegree(i);
  }

  [[nodiscard]] std::shared_ptr<baf::KnotVector> GetKnotVector(int i) const {
    return parameter_space_->GetKnotVector(i);
  }

  [[nodiscard]] double GetKnotVectorRange(int direction) const {
    return parameter_space_->GetKnotVectorRange(direction);
  }

  [[nodiscard]] int GetNumberOfControlPoints() const {
    return GetPhysicalSpace()->GetNumberOfControlPoints();
  }

  std::array<int, DIM> GetPointsPerDirection() const {
    return GetPhysicalSpace()->GetPointsPerDirection();
  }

  [[nodiscard]] int GetPointDim() const {
    return GetPhysicalSpace()->GetDimension();
  }

  double GetControlPoint(std::array<int, DIM> indices, int dimension) const {
    return GetPhysicalSpace()->GetControlPoint(indices).GetValue(dimension);
  }

  baf::ControlPoint GetControlPoint(std::array<int, DIM> indices) const {
    return GetPhysicalSpace()->GetControlPoint(indices);
  }

  [[nodiscard]] double GetExpansion() const {
    return GetPhysicalSpace()->GetExpansion();
  }

  double GetWeight(std::array<int, DIM> indices) const {
    return GetPhysicalSpace()->GetWeight(indices);
  }

  void InsertKnot(ParamCoord knot, int dimension, size_t multiplicity = 1) {
    KnotSpan knot_span = GetKnotVector(dimension)->GetKnotSpan(knot);
    Degree degree = GetDegree(dimension);
    for (size_t i = 1; i <= multiplicity; ++i) {
      if (GetKnotVector(dimension)->IsLastKnot(knot)) {
        knot_span = knot_span + KnotSpan{1};
      }
      auto last = knot_span.get() - GetKnotVector(dimension)->GetMultiplicity(knot);
      auto first = knot_span.get() - degree.get() + i;
      std::vector<double> scaling;
      for (auto j = static_cast<size_t>(first); j <= last; ++j) {
        ParamCoord low_knot = GetKnotVector(dimension)->GetKnot(j);
        ParamCoord upper_knot = GetKnotVector(dimension)->GetKnot(j + degree.get() - i + 1);
        scaling.emplace_back((knot.get() - low_knot.get()) / (upper_knot.get() - low_knot.get()));
      }
      this->AdjustControlPoints(scaling, static_cast<int>(first), static_cast<int>(last), dimension);
    }
    for (size_t i = 0; i < multiplicity; ++i) {
      parameter_space_->InsertKnot(knot, dimension);
    }
  }

  void RefineKnots(const std::vector<ParamCoord> &new_knots, int dimension) {
    for (const auto &knot : new_knots) {
      this->InsertKnot(knot, dimension);
    }
  }

  size_t RemoveKnot(ParamCoord knot, int dimension, double tolerance, size_t multiplicity = 1) {
    size_t count = 0;
    for (; count < multiplicity; ++count) {
      KnotSpan knot_span = GetKnotVector(dimension)->GetKnotSpan(knot);
      Degree degree = GetDegree(dimension);
      auto last = knot_span.get() - GetKnotVector(dimension)->GetMultiplicity(knot);
      auto first = knot_span.get() - degree.get();
      std::vector<double> scaling;
      for (auto i = static_cast<size_t>(first); i <= last; ++i) {
        ParamCoord low_knot = GetKnotVector(dimension)->GetKnot(i);
        ParamCoord upper_knot = GetKnotVector(dimension)->GetKnot(i + degree.get() + 1);
        scaling.emplace_back((knot.get() - low_knot.get()) / (upper_knot.get() - low_knot.get()));
      }
      bool is_removed = this->RemoveControlPoints(scaling, first, static_cast<int>(last), dimension, tolerance);
      if (!is_removed) break;
      parameter_space_->RemoveKnot(knot, dimension);
    }
    return count;
  }

  void ElevateDegreeForDimension(int dimension) {
    uint64_t num_bezier_segments = GetKnotVector(dimension)->GetNumberOfDifferentKnots() - 1;
    std::vector<double> alpha = ComputeBezierDegreeElevationCoeffients(dimension);
    auto diff = ProduceBezierSegments(dimension);
    std::vector<std::vector<baf::ControlPoint>> bezier_segments(num_bezier_segments);
    for (uint64_t i = 0; i < num_bezier_segments; ++i) {
      bezier_segments[i] = DegreeElevateBezierSegment(GetBezierSegment(dimension, i), alpha, dimension);
    }
    std::array<int, DIM> cps_per_dir = GetPointsPerDirection();
    cps_per_dir[dimension] += (GetKnotVector(dimension)->GetNumberOfDifferentKnots() - 1);
    util::MultiIndexHandler<DIM> cp_handler(cps_per_dir);
    GetPhysicalSpace()->SetNumberOfPoints(dimension, cps_per_dir[dimension]);
    int delta_num_cps = cp_handler.Get1DLength() - GetPhysicalSpace()->GetNumberOfControlPoints();
    GetPhysicalSpace()->AddControlPoints(delta_num_cps);
    parameter_space_->ElevateDegree(dimension);
    parameter_space_->IncrementMultiplicityOfAllKnots(dimension);
    SetNewBezierSegmentControlPoints(bezier_segments, dimension);
    RemoveBezierKnots(diff, dimension);
  }

  bool ReduceDegreeForDimension(int dimension, double tolerance) {
    std::vector<int> diff = ProduceBezierSegments(dimension);
    uint64_t num_bezier_segments = GetKnotVector(dimension)->GetNumberOfDifferentKnots() - 1;
    std::vector<std::vector<baf::ControlPoint>> bezier_segments;
    for (uint64_t i = 0; i < num_bezier_segments; ++i) {
      bool successful;
      bezier_segments.emplace_back(
          DegreeReduceBezierSegment(GetBezierSegment(dimension, i), tolerance, dimension, &successful));
      if (!successful) {
        RemoveBezierKnots(diff, dimension);
        return false;
      }
    }
    std::array<int, DIM> cps_per_dir = GetPointsPerDirection();
    cps_per_dir[dimension] -= num_bezier_segments;
    GetPhysicalSpace()->SetNumberOfPoints(dimension, cps_per_dir[dimension]);
    util::MultiIndexHandler<DIM> point_handler(cps_per_dir);
    int delta_num_cps = GetPhysicalSpace()->GetNumberOfControlPoints() - point_handler.Get1DLength();
    GetPhysicalSpace()->RemoveControlPoints(delta_num_cps);
    parameter_space_->ReduceDegree(dimension);
    parameter_space_->DecrementMultiplicityOfAllKnots(dimension);
    SetNewBezierSegmentControlPoints(bezier_segments, dimension);
    RemoveBezierKnots(diff, dimension);
    return true;
  }

  virtual void AdjustControlPoints(std::vector<double> scaling, int first, int last, int dimension) = 0;
  virtual bool RemoveControlPoints(std::vector<double> scaling,
                                   int first, int last, int dimension, double tolerance) = 0;

 protected:
  void ThrowIfParametricCoordinateOutsideKnotVectorRange(std::array<ParamCoord, DIM> param_coord) const {
    parameter_space_->ThrowIfParametricCoordinateOutsideKnotVectorRange(param_coord);
  }

  virtual double GetEvaluatedControlPoint(std::array<ParamCoord, DIM> param_coord,
                                          std::array<int, DIM> indices,
                                          int dimension) const = 0;

  virtual double GetEvaluatedDerivativeControlPoint(std::array<ParamCoord, DIM> param_coord,
                                                    std::array<int, DIM> derivative,
                                                    std::array<int, DIM> indices,
                                                    int dimension) const = 0;

  virtual std::shared_ptr<spl::PhysicalSpace<DIM>> GetPhysicalSpace() const = 0;

  std::array<int, DIM> GetArrayOfFirstNonZeroBasisFunctions(std::array<ParamCoord, DIM> param_coord) const {
    return parameter_space_->GetArrayOfFirstNonZeroBasisFunctions(param_coord);
  }

  std::array<int, DIM> GetNumberOfBasisFunctionsToEvaluate() const {
    std::array<int, DIM> total_length;
    for (int i = 0; i < DIM; ++i) {
      total_length[i] = parameter_space_->GetDegree(i).get() + 1;
    }
    return total_length;
  }

  [[nodiscard]] std::vector<double> ComputeBezierDegreeElevationCoeffients(int dimension) const {
    std::vector<double> alpha;
    for (int i = 0; i < GetDegree(dimension).get() + 2; ++i) {
      alpha.push_back(static_cast<double>(i) / (GetDegree(dimension).get() + 1));
    }
    return alpha;
  }

  std::vector<int> ProduceBezierSegments(int dimension) {
    size_t position = GetKnotVector(dimension)->GetNumberOfKnots() - GetDegree(dimension).get() - 2;
    std::vector<int> diff(GetKnotVector(dimension)->GetNumberOfDifferentKnots() - 2, GetDegree(dimension).get() - 1);

    for (size_t i = GetKnotVector(dimension)->GetNumberOfDifferentKnots() - 2; i > 0; --i, --position) {
      while (GetKnotVector(dimension)->GetKnot(position) == GetKnotVector(dimension)->GetKnot(position - 1)) {
        --position;
        --diff[i - 1];
      }
      InsertKnot(GetKnotVector(dimension)->GetKnot(position), dimension, static_cast<size_t>(diff[i - 1]));
    }
    return diff;
  }

  void RemoveBezierKnots(std::vector<int> diff, int dimension) {
    size_t position = GetKnotVector(dimension)->GetNumberOfKnots() - GetDegree(dimension).get() - 2;
    for (size_t i = GetKnotVector(dimension)->GetNumberOfDifferentKnots() - 2; i > 0; --i, --position) {
      while (GetKnotVector(dimension)->GetKnot(position) == GetKnotVector(dimension)->GetKnot(position - 1)) {
        --position;
      }
      RemoveKnot(GetKnotVector(dimension)->GetKnot(position), dimension, util::NumericSettings<double>::kEpsilon(),
                 static_cast<size_t>(diff[i - 1]));
    }
  }

  [[nodiscard]] std::vector<baf::ControlPoint> GetBezierSegment(int dimension, int segment) const {
    util::MultiIndexHandler<DIM> point_handler(GetPointsPerDirection());
    int width = GetDegree(dimension).get() + 1;
    int segment_length = GetNumberOfControlPoints() / GetPointsPerDirection()[dimension];
    std::vector<baf::ControlPoint> bezier_cps(static_cast<size_t>(width * segment_length),
                                              baf::ControlPoint(GetPointDim() + 1));
    for (int i = 0; i < point_handler.Get1DLength(); ++i, ++point_handler) {
      auto index_in_dir = point_handler[dimension];
      if (index_in_dir >= segment * width - segment && index_in_dir < (segment + 1) * width - segment) {
        auto index = ((index_in_dir + segment) % width + point_handler.ExtractDimension(dimension) * width);
        double weight = GetWeight(point_handler.GetIndices());
        for (int j = 0; j < GetPointDim(); ++j) {
          bezier_cps[index].SetValue(j, GetControlPoint(point_handler.GetIndices()).GetValue(j) * weight);
        }
        bezier_cps[index].SetValue(GetPointDim(), weight);
      }
    }
    return bezier_cps;
  }

  [[nodiscard]] std::vector<baf::ControlPoint> DegreeElevateBezierSegment(
      const std::vector<baf::ControlPoint> &bezier_cps, std::vector<double> alpha, int dimension) const {
    int width = GetDegree(dimension).get() + 2;
    int segment_length = GetNumberOfControlPoints() / GetPointsPerDirection()[dimension];
    util::MultiIndexHandler<2> point_handler({width, segment_length});
    std::vector<baf::ControlPoint> new_cps;
    for (int k = 0; k < point_handler.Get1DLength(); ++k, ++point_handler) {
      auto index = point_handler[0];
      auto pos = point_handler.Get1DIndex() - point_handler[1];
      if (alpha[index] == 0) {
        new_cps.emplace_back(bezier_cps[pos] * (1 - alpha[index]));
      } else if (alpha[index] == 1) {
        new_cps.emplace_back(bezier_cps[pos - 1] * alpha[index]);
      } else {
        new_cps.emplace_back(bezier_cps[pos] * (1 - alpha[index]) + bezier_cps[pos - 1] * alpha[index]);
      }
    }
    return new_cps;
  }

  std::vector<baf::ControlPoint> DegreeReduceBezierSegment(const std::vector<baf::ControlPoint> &bezier_cps,
      double tolerance, int dimension, bool* successful) {
    double max_error_bound = 0.0;
    int new_degree_in_dimension = GetDegree(dimension).get() - 1;
    int width = GetDegree(dimension).get() + 1;
    std::array<int, DIM> points_per_dir = GetPointsPerDirection();
    points_per_dir[dimension] = width;
    util::MultiIndexHandler<DIM> cp_handler_old(points_per_dir);
    points_per_dir[dimension]--;
    util::MultiIndexHandler<DIM> cp_handler_new(points_per_dir);
    std::vector<baf::ControlPoint> control_points_new(cp_handler_new.Get1DLength(), baf::ControlPoint({0.0}));
    auto indices_new = cp_handler_new.Get1DIndicesForFixedDimension(dimension, 0);
    auto indices_old = cp_handler_old.Get1DIndicesForFixedDimension(dimension, 0);
    for (uint64_t i = 0; i < indices_new.size(); ++i) {
      control_points_new[indices_new[i]] = bezier_cps[indices_old[i]];
    }
    indices_new = cp_handler_new.Get1DIndicesForFixedDimension(dimension, new_degree_in_dimension);
    indices_old = cp_handler_old.Get1DIndicesForFixedDimension(dimension, new_degree_in_dimension + 1);
    for (uint64_t i = 0; i < indices_new.size(); ++i) {
      control_points_new[indices_new[i]] = bezier_cps[indices_old[i]];
    }
    int r = new_degree_in_dimension / 2;
    int end_of_first_loop = r;
    if ((new_degree_in_dimension + 1) % 2 != 0) end_of_first_loop = r - 1;
    for (int i = 1; i <= end_of_first_loop; ++i) {
      auto indices_new_minus_one = cp_handler_new.Get1DIndicesForFixedDimension(dimension, i - 1);
      indices_new = cp_handler_new.Get1DIndicesForFixedDimension(dimension, i);
      indices_old = cp_handler_old.Get1DIndicesForFixedDimension(dimension, i);
      double alpha = static_cast<double>(i) / static_cast<double>(new_degree_in_dimension + 1);
      for (uint64_t j = 0; j < indices_new.size(); ++j) {
        control_points_new[indices_new[j]] = (bezier_cps[indices_old[j]]
            - control_points_new[indices_new_minus_one[j]] * alpha) * (1 / (1 - alpha));
      }
    }
    int lower = static_cast<int>(r + 1);
    int upper = new_degree_in_dimension - 1;
    for (int i = upper; i >= lower; --i) {
      auto indices_new_plus_one = cp_handler_new.Get1DIndicesForFixedDimension(dimension, i + 1);
      indices_new = cp_handler_new.Get1DIndicesForFixedDimension(dimension, i);
      auto indices_old_plus_one = cp_handler_old.Get1DIndicesForFixedDimension(dimension, i + 1);
      double alpha = static_cast<double>(i + 1) / static_cast<double>(new_degree_in_dimension + 1);
      for (uint64_t j = 0; j < indices_new.size(); ++j) {
        control_points_new[indices_new[j]] = (bezier_cps[indices_old_plus_one[j]]
            - control_points_new[indices_new_plus_one[j]] * (1 - alpha)) * (1 / alpha);
      }
    }
    if ((new_degree_in_dimension + 1) % 2 != 0) {
      auto indices_new_r = cp_handler_new.Get1DIndicesForFixedDimension(dimension, r);
      auto indices_new_r_minus_one = cp_handler_new.Get1DIndicesForFixedDimension(dimension, r - 1);
      auto indices_new_r_plus_one = cp_handler_new.Get1DIndicesForFixedDimension(dimension, r + 1);
      auto indices_old_r = cp_handler_old.Get1DIndicesForFixedDimension(dimension, r);
      auto indices_old_r_plus_one = cp_handler_old.Get1DIndicesForFixedDimension(dimension, r + 1);
      double alpha_r = static_cast<double>(r) / static_cast<double>((new_degree_in_dimension + 1));
      double alpha_r_plus_one = static_cast<double>(r + 1) / static_cast<double>((new_degree_in_dimension + 1));
      for (uint64_t j = 0; j < indices_new.size(); ++j) {
        baf::ControlPoint P_r_L = (bezier_cps[indices_old_r[j]]
            - control_points_new[indices_new_r_minus_one[j]] * alpha_r) * (1 / (1 - alpha_r));
        baf::ControlPoint P_r_R = (bezier_cps[indices_old_r_plus_one[j]]
            - control_points_new[indices_new_r_plus_one[j]] * (1 - alpha_r_plus_one)) * (1 / alpha_r_plus_one);
        control_points_new[indices_new_r[j]] = (P_r_L + P_r_R) * 0.5;
        if ((P_r_L - P_r_R).GetEuclideanNorm() > max_error_bound) max_error_bound = (P_r_L - P_r_R).GetEuclideanNorm();
      }
    }
    if ((new_degree_in_dimension + 1) % 2 == 0) {
      auto indices_new_r = cp_handler_new.Get1DIndicesForFixedDimension(dimension, r);
      auto indices_new_r_plus_one = cp_handler_new.Get1DIndicesForFixedDimension(dimension, r + 1);
      auto indices_old_r_plus_one = cp_handler_old.Get1DIndicesForFixedDimension(dimension, r + 1);
      for (uint64_t j = 0; j < indices_new.size(); ++j) {
        double current_max_error_bound =
            (bezier_cps[indices_old_r_plus_one[j]] - ((control_points_new[indices_new_r[j]]
                + control_points_new[indices_new_r_plus_one[j]]) * 0.5)).GetEuclideanNorm();
        if (current_max_error_bound > max_error_bound) max_error_bound = current_max_error_bound;
      }
    }
    *successful = max_error_bound < tolerance;
    return control_points_new;
  }

  virtual void SetNewBezierPoint(baf::ControlPoint control_point, double weight, std::array<int, DIM> indices) = 0;

  void SetNewBezierSegmentControlPoints(const std::vector<std::vector<baf::ControlPoint>> &bezier_segments,
                                        int dimension) {
    util::MultiIndexHandler<DIM> point_handler(GetPointsPerDirection());
    int width = GetDegree(dimension).get() + 1;
    for (int i = 0; i < point_handler.Get1DLength(); ++i, ++point_handler) {
      int index_in_dir = point_handler[dimension];
      int segment = index_in_dir / (width - 1);
      int index = index_in_dir % (width - 1) + point_handler.ExtractDimension(dimension) * width;
      if (static_cast<size_t>(segment) > bezier_segments.size() - 1) {
        segment--;
        index = (width - 1) + point_handler.ExtractDimension(dimension) * width;
      }
      baf::ControlPoint cps(GetPointDim());
      double weight = bezier_segments[segment][index].GetValue(GetPointDim());
      for (int j = 0; j < GetPointDim(); ++j) {
        cps.SetValue(j, bezier_segments[segment][index].GetValue(j) / weight);
      }
      this->SetNewBezierPoint(cps, weight, point_handler.GetIndices());
    }
  }

  std::shared_ptr<ParameterSpace<DIM>> parameter_space_;
};
}  // namespace spl

#endif  // SRC_SPL_SPLINE_H_
