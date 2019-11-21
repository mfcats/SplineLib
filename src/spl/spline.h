/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#ifndef SRC_SPL_SPLINE_H_
#define SRC_SPL_SPLINE_H_

#include <algorithm>
#include <array>
#include <functional>
#include <numeric>
#include <sstream>
#include <utility>
#include <vector>

#include "src/spl/control_point.h"
#include "src/baf/knot_vector.h"
#include "src/util/multi_index_handler.h"
#include "src/spl/parameter_space.h"
#include "src/spl/physical_space.h"

namespace splinelib::src::spl {
template<int PARAMETRIC_DIMENSIONALITY>
class Spline {
 public:
  virtual ~Spline() = default;
  Spline() = default;
  Spline(baf::KnotVectors<PARAMETRIC_DIMENSIONALITY> knot_vector,
      std::array<Degree, PARAMETRIC_DIMENSIONALITY> degree) {
    parameter_space_ = std::make_shared<ParameterSpace<PARAMETRIC_DIMENSIONALITY>>(
        ParameterSpace<PARAMETRIC_DIMENSIONALITY>(knot_vector, degree));
  }
  explicit Spline(std::shared_ptr<ParameterSpace<PARAMETRIC_DIMENSIONALITY>> parameter_space) {
    parameter_space_ = parameter_space;
  }
  Spline(const Spline<PARAMETRIC_DIMENSIONALITY> &spline) {
    ParameterSpace<PARAMETRIC_DIMENSIONALITY> parameter_space(*spline.parameter_space_);
    parameter_space_ = std::make_shared<ParameterSpace<PARAMETRIC_DIMENSIONALITY>>(parameter_space);
  }

  Spline(Spline<PARAMETRIC_DIMENSIONALITY> &&other) = delete;
  Spline & operator=(const Spline<PARAMETRIC_DIMENSIONALITY> &rhs) = delete;
  Spline & operator=(Spline<PARAMETRIC_DIMENSIONALITY> &&rhs) = delete;

  virtual std::vector<double> Evaluate(std::array<ParametricCoordinate, PARAMETRIC_DIMENSIONALITY> param_coord,
                                       const std::vector<int> &dimensions) const {
    this->ThrowIfParametricCoordinateOutsideKnotVectorRange(param_coord);
    auto first_non_zero = GetArrayOfFirstNonZeroBasisFunctions(param_coord);
    util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> basisFunctionHandler(
        this->GetNumberOfBasisFunctionsToEvaluate());
    std::vector<double> evaluated_point(dimensions.size(), 0);
    for (int i = 0; i < basisFunctionHandler.GetNumberOfTotalMultiIndices(); ++i, basisFunctionHandler++) {
      auto indices = basisFunctionHandler.GetCurrentIndex();
      std::transform(indices.begin(), indices.end(), first_non_zero.begin(), indices.begin(), std::plus<>());
      for (uint64_t j = 0; j < dimensions.size(); ++j) {
        evaluated_point[j] += GetEvaluatedControlPoint(param_coord, indices, dimensions[j]);
      }
    }
    return evaluated_point;
  }

  virtual std::vector<double> EvaluateDerivative(
      std::array<ParametricCoordinate, PARAMETRIC_DIMENSIONALITY> param_coord, const std::vector<int> &dimensions,
      std::array<int, PARAMETRIC_DIMENSIONALITY> derivative) const {
    this->ThrowIfParametricCoordinateOutsideKnotVectorRange(param_coord);

    auto first_non_zero = this->GetArrayOfFirstNonZeroBasisFunctions(param_coord);
    util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> basisFunctionHandler(
        this->GetNumberOfBasisFunctionsToEvaluate());
    std::vector<double> evaluated_point(dimensions.size(), 0);

    for (int i = 0; i < basisFunctionHandler.GetNumberOfTotalMultiIndices(); ++i, basisFunctionHandler++) {
      auto indices = basisFunctionHandler.GetCurrentIndex();
      std::transform(indices.begin(), indices.end(), first_non_zero.begin(), indices.begin(), std::plus<>());
      for (uint64_t j = 0; j < dimensions.size(); ++j) {
        evaluated_point[j] += GetEvaluatedDerivativeControlPoint(param_coord, derivative, indices, dimensions[j]);
      }
    }
    return evaluated_point;
  }

  std::vector<double> EvaluateAllNonZeroBasisFunctions(int direction, ParametricCoordinate param_coord) const {
    return parameter_space_->EvaluateAllNonZeroBasisFunctions(direction, param_coord);
  }

  std::vector<double> EvaluateAllNonZeroBasisFunctionDerivatives(int direction,
                                                                 ParametricCoordinate param_coord,
                                                                 int derivative) const {
    return parameter_space_->EvaluateAllNonZeroBasisFunctionDerivatives(direction, param_coord, derivative);
  }

  bool AreGeometricallyEqual(const spl::Spline<PARAMETRIC_DIMENSIONALITY> &rhs,
                             double tolerance = util::numeric_settings::GetEpsilon<double>()) const {
    double number = ceil(pow(100, 1.0 / PARAMETRIC_DIMENSIONALITY));
    std::array<int, PARAMETRIC_DIMENSIONALITY> points = std::array<int, PARAMETRIC_DIMENSIONALITY>();
    for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; ++i) {
      points[i] = number + 1;
    }
    std::vector<int> dimensions;
    dimensions.reserve(GetPointDim());
    for (int i = 0; i < GetPointDim(); ++i) {
      dimensions.emplace_back(i);
    }
    util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> point_handler(points);
    for (int i = 0; i < point_handler.GetNumberOfTotalMultiIndices(); ++i, ++point_handler) {
      std::array<ParametricCoordinate, PARAMETRIC_DIMENSIONALITY> param_coord{};
      for (int dim = 0; dim < PARAMETRIC_DIMENSIONALITY; ++dim) {
        double span = GetKnotVector(dim)->GetLastKnot().Get() - GetKnotVector(dim)->GetKnot(0).Get();
        param_coord[dim] = ParametricCoordinate(span / number * point_handler[Dimension{dim}]);
      }
      std::vector<double> evaluate_this = Evaluate(param_coord, dimensions);
      std::vector<double> evaluate_rhs = rhs.Evaluate(param_coord, dimensions);
      if (util::vector_utils::ComputeDistance(evaluate_this, evaluate_rhs) > tolerance) {
        return false;
      }
    }
    return true;
  }

  Degree GetDegree(int i) const {
    return parameter_space_->GetDegree(i);
  }

  std::shared_ptr<baf::KnotVector> GetKnotVector(int i) const {
    return parameter_space_->GetKnotVector(i);
  }

  double GetKnotVectorRange(int direction) const {
    return parameter_space_->GetKnotVectorRange(direction);
  }

  int GetNumberOfControlPoints() const {
    return GetPhysicalSpace()->GetNumberOfControlPoints();
  }

  std::array<int, PARAMETRIC_DIMENSIONALITY> GetPointsPerDirection() const {
    return GetPhysicalSpace()->GetPointsPerDirection();
  }

  int GetPointDim() const {
    return GetPhysicalSpace()->GetDimension();
  }

  double GetControlPoint(std::array<int, PARAMETRIC_DIMENSIONALITY> indices, int dimension) const {
    return GetPhysicalSpace()->GetControlPoint(indices).GetValueForDimension(Dimension{dimension});
  }

  double GetControlPoint(int index_1d, int dimension) const {
    return GetPhysicalSpace()->GetControlPoint(index_1d).GetValueForDimension(Dimension{dimension});
  }

  baf::ControlPoint GetControlPoint(std::array<int, PARAMETRIC_DIMENSIONALITY> indices) const {
    return GetPhysicalSpace()->GetControlPoint(indices);
  }

  baf::ControlPoint GetControlPoint(int index_1d) const {
    return GetPhysicalSpace()->GetControlPoint(index_1d);
  }

  double GetExpansion() const {
    return GetPhysicalSpace()->GetExpansion();
  }

  double GetWeight(std::array<int, PARAMETRIC_DIMENSIONALITY> indices) const {
    return GetPhysicalSpace()->GetWeight(indices);
  }

  double GetWeight(int index_1d) const {
    return GetPhysicalSpace()->GetWeight(index_1d);
  }

  void InsertKnot(ParametricCoordinate knot, int dimension, size_t multiplicity = 1) {
    KnotSpan knot_span = GetKnotVector(dimension)->GetKnotSpan(knot);
    Degree degree = GetDegree(dimension);
    for (size_t i = 1; i <= multiplicity; ++i) {
      if (GetKnotVector(dimension)->IsLastKnot(knot)) knot_span = knot_span + KnotSpan{1};
      auto last = knot_span.Get() - GetKnotVector(dimension)->GetMultiplicity(knot).Get();
      auto first = knot_span.Get() - degree.Get() + i;
      std::vector<double> scaling;
      for (int j = first; j <= last; ++j) {
        ParametricCoordinate low_knot = GetKnotVector(dimension)->GetKnot(j);
        ParametricCoordinate upper_knot = GetKnotVector(dimension)->GetKnot(j + degree.Get() - i + 1);
        scaling.emplace_back((knot.Get() - low_knot.Get()) / (upper_knot.Get() - low_knot.Get()));
      }
      this->AdjustControlPoints(scaling, static_cast<int>(first), static_cast<int>(last), dimension);
    }
    for (size_t i = 0; i < multiplicity; ++i) {
      parameter_space_->InsertKnot(knot, dimension);
    }
  }

  void RefineKnots(std::vector<ParametricCoordinate> new_knots, int dimension) {
    for (const auto &knot : new_knots) {
      this->InsertKnot(knot, dimension);
    }
  }

  size_t RemoveKnot(ParametricCoordinate knot, int dimension, double tolerance, size_t multiplicity = 1) {
    size_t count = 0;
    for (; count < multiplicity; ++count) {
      KnotSpan knot_span = GetKnotVector(dimension)->GetKnotSpan(knot);
      Degree degree = GetDegree(dimension);
      auto last = knot_span.Get() - GetKnotVector(dimension)->GetMultiplicity(knot).Get();
      auto first = knot_span.Get() - degree.Get();
      std::vector<double> scaling;
      for (int i = first; i <= last; ++i) {
        ParametricCoordinate low_knot = GetKnotVector(dimension)->GetKnot(i);
        ParametricCoordinate upper_knot = GetKnotVector(dimension)->GetKnot(i + degree.Get() + 1);
        scaling.emplace_back((knot.Get() - low_knot.Get()) / (upper_knot.Get() - low_knot.Get()));
      }
      bool is_removed = this->RemoveControlPoints(scaling, first, static_cast<int>(last), dimension, tolerance);
      if (!is_removed) break;
      parameter_space_->RemoveKnot(knot, dimension);
    }
    return count;
  }

  void ElevateDegreeForDimension(int dimension) {
    uint64_t num_bezier_segments = GetKnotVector(dimension)->GetNumberOfDifferentKnots() - 1;
    std::vector<double> alpha = ComputeBezierDegreeElevationCoefficients(dimension);
    auto diff = ProduceBezierSegments(dimension);
    std::vector<std::vector<baf::ControlPoint>> bezier_segments(num_bezier_segments);
    for (uint64_t i = 0; i < num_bezier_segments; ++i) {
      bezier_segments[i] = DegreeElevateBezierSegment(GetBezierSegment(dimension, i), alpha, dimension);
    }
    std::array<int, PARAMETRIC_DIMENSIONALITY> cps_per_dir = GetPointsPerDirection();
    cps_per_dir[dimension] += (GetKnotVector(dimension)->GetNumberOfDifferentKnots() - 1);
    util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> cp_handler(cps_per_dir);
    GetPhysicalSpace()->SetNumberOfPoints(dimension, cps_per_dir[dimension]);
    int delta_num_cps = cp_handler.GetNumberOfTotalMultiIndices() - GetPhysicalSpace()->GetNumberOfControlPoints();
    GetPhysicalSpace()->AddControlPoints(delta_num_cps);
    parameter_space_->ElevateDegree(dimension);
    parameter_space_->IncrementMultiplicityOfAllKnots(dimension);
    SetNewBezierSegmentControlPoints(bezier_segments, dimension);
    RemoveBezierKnots(diff, dimension);
  }

  bool ReduceDegreeForDimension(int dimension, double tolerance = util::numeric_settings::GetEpsilon<double>()) {
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
    std::array<int, PARAMETRIC_DIMENSIONALITY> cps_per_dir = GetPointsPerDirection();
    cps_per_dir[dimension] -= num_bezier_segments;
    GetPhysicalSpace()->SetNumberOfPoints(dimension, cps_per_dir[dimension]);
    util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> point_handler(cps_per_dir);
    int delta_num_cps = GetPhysicalSpace()->GetNumberOfControlPoints() - point_handler.GetNumberOfTotalMultiIndices();
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
  void ThrowIfParametricCoordinateOutsideKnotVectorRange(
      std::array<ParametricCoordinate, PARAMETRIC_DIMENSIONALITY> param_coord) const {
    parameter_space_->ThrowIfParametricCoordinateOutsideKnotVectorRange(param_coord);
  }

  virtual double GetEvaluatedControlPoint(std::array<ParametricCoordinate, PARAMETRIC_DIMENSIONALITY> param_coord,
                                          std::array<int, PARAMETRIC_DIMENSIONALITY> indices,
                                          int dimension) const = 0;

  virtual double GetEvaluatedDerivativeControlPoint(
      std::array<ParametricCoordinate, PARAMETRIC_DIMENSIONALITY> param_coord,
      std::array<int, PARAMETRIC_DIMENSIONALITY> derivative, std::array<int, PARAMETRIC_DIMENSIONALITY> indices,
      int dimension) const = 0;

  virtual std::shared_ptr<spl::PhysicalSpace<PARAMETRIC_DIMENSIONALITY>> GetPhysicalSpace() const = 0;

  std::array<int, PARAMETRIC_DIMENSIONALITY> GetArrayOfFirstNonZeroBasisFunctions(
      std::array<ParametricCoordinate, PARAMETRIC_DIMENSIONALITY> param_coord) const {
    return parameter_space_->GetArrayOfFirstNonZeroBasisFunctions(param_coord);
  }

  std::array<int, PARAMETRIC_DIMENSIONALITY> GetNumberOfBasisFunctionsToEvaluate() const {
    std::array<int, PARAMETRIC_DIMENSIONALITY> total_length{};
    for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; ++i) {
      total_length[i] = parameter_space_->GetDegree(i).Get() + 1;
    }
    return total_length;
  }

  std::vector<double> ComputeBezierDegreeElevationCoefficients(int dimension) const {
    std::vector<double> alpha;
    for (int current_bezier_point = 0; current_bezier_point < GetDegree(dimension).Get() + 2; ++current_bezier_point) {
      alpha.emplace_back(current_bezier_point / (GetDegree(dimension).Get() + 1.0));
    }
    return alpha;
  }

  std::vector<int> ProduceBezierSegments(int dimension) {
    uint64_t current_knot = GetDegree(dimension).Get() + 1;
    std::vector<int> diff(GetKnotVector(dimension)->GetNumberOfDifferentKnots() - 2, GetDegree(dimension).Get() - 1);
    for (int current_knot_span = 0; current_knot_span < GetKnotVector(dimension)->GetNumberOfDifferentKnots() - 2;
         ++current_knot_span, ++current_knot) {
      while (GetKnotVector(dimension)->GetKnot(current_knot) == GetKnotVector(dimension)->GetKnot(current_knot + 1)) {
        ++current_knot;
        --diff[current_knot_span];
      }
      InsertKnot(GetKnotVector(dimension)->GetKnot(current_knot), dimension,
                 static_cast<size_t>(diff[current_knot_span]));
      current_knot += diff[current_knot_span];
    }
    return diff;
  }

  void RemoveBezierKnots(std::vector<int> diff, int dimension) {
    auto current_knot = GetDegree(dimension).Get() + 1;
    for (int current_knot_span = 0; current_knot_span < GetKnotVector(dimension)->GetNumberOfDifferentKnots() - 2;
         ++current_knot_span, ++current_knot) {
      RemoveKnot(GetKnotVector(dimension)->GetKnot(current_knot), dimension,
                 util::numeric_settings::GetEpsilon<double>(), static_cast<size_t>(diff[current_knot_span]));
      while (GetKnotVector(dimension)->GetKnot(current_knot) == GetKnotVector(dimension)->GetKnot(current_knot + 1)) {
        ++current_knot;
      }
    }
  }

  std::vector<baf::ControlPoint> GetBezierSegment(int dimension, int segment) const {
    util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> point_handler(GetPointsPerDirection());
    int width = GetDegree(dimension).Get() + 1;
    int segment_length = GetNumberOfControlPoints() / GetPointsPerDirection()[dimension];
    std::vector<baf::ControlPoint> bezier_cps(static_cast<size_t>(width * segment_length),
                                              baf::ControlPoint(GetPointDim() + 1));
    for (int i = 0; i < point_handler.GetNumberOfTotalMultiIndices(); ++i, ++point_handler) {
      auto index_in_dir = point_handler[Dimension{dimension}];
      if (index_in_dir >= segment * width - segment && index_in_dir < (segment + 1) * width - segment) {
        auto index = ((index_in_dir + segment) % width + point_handler.CollapseDimension(Dimension{dimension}) * width);
        double weight = GetWeight(point_handler.GetCurrentIndex());
        for (int j = 0; j < GetPointDim(); ++j) {
          bezier_cps[index].SetValue(j, GetControlPoint(
              point_handler.GetCurrentIndex()).GetValueForDimension(Dimension{j}) * weight);
        }
        bezier_cps[index].SetValue(GetPointDim(), weight);
      }
    }
    return bezier_cps;
  }

  std::vector<baf::ControlPoint> DegreeElevateBezierSegment(const std::vector<baf::ControlPoint> &bezier_cps,
                                                            std::vector<double> alpha, int dimension) const {
    int width = GetDegree(dimension).Get() + 2;
    int segment_length = GetNumberOfControlPoints() / GetPointsPerDirection()[dimension];
    util::MultiIndexHandler<2> cp_handler_new({width, segment_length});
    util::MultiIndexHandler<2> cp_handler_old({width - 1, segment_length});
    std::vector<baf::ControlPoint> new_cps(cp_handler_new.GetNumberOfTotalMultiIndices(),
                                           baf::ControlPoint(GetPointDim() + 1));
    for (int i = 0; i < segment_length; ++i) {
      new_cps[cp_handler_new.Get1DIndex({0, i})] = bezier_cps[cp_handler_old.Get1DIndex({0, i})];
      new_cps[cp_handler_new.Get1DIndex({width - 1, i})] = bezier_cps[cp_handler_old.Get1DIndex({width - 2, i})];
      for (int j = 1; j < GetDegree(dimension).Get(); ++j) {
        new_cps[cp_handler_new.Get1DIndex({j, i})] = bezier_cps[cp_handler_old.Get1DIndex({j, i})] * (1 - alpha[j]) +
            alpha[j] * bezier_cps[cp_handler_old.Get1DIndex({j - 1, i})] /*  * alpha[j]*/;
      }
    }
    return new_cps;
  }

  std::vector<baf::ControlPoint> DegreeReduceBezierSegment(const std::vector<baf::ControlPoint> &bezier_cps,
                                                           double tolerance, int dimension, bool* successful) {
    int segment_length = GetNumberOfControlPoints() / GetPointsPerDirection()[dimension];
    util::MultiIndexHandler<2> cp_handler_old({GetDegree(dimension).Get() + 1, segment_length});
    util::MultiIndexHandler<2> cp_handler_new({GetDegree(dimension).Get(), segment_length});
    double max_error_bound = 0.0;
    int new_degree = GetDegree(dimension).Get() - 1;
    std::vector<baf::ControlPoint> new_cps(cp_handler_new.GetNumberOfTotalMultiIndices(),
                                           baf::ControlPoint(GetPointDim() + 1));
    int r = new_degree / 2;
    double alpha_r = static_cast<double>(r) / static_cast<double>((new_degree + 1));
    double alpha_r_plus = static_cast<double>(r + 1) / static_cast<double>((new_degree + 1));
    int end_of_first_loop = r;
    if ((new_degree + 1) % 2 != 0) end_of_first_loop = r - 1;
    for (int i = 0; i < segment_length; ++i) {
      new_cps[cp_handler_new.Get1DIndex({0, i})] = bezier_cps[cp_handler_old.Get1DIndex({0, i})];
      new_cps[cp_handler_new.Get1DIndex({new_degree, i})] = bezier_cps[cp_handler_old.Get1DIndex({new_degree + 1, i})];
      for (int j = 1; j <= end_of_first_loop; ++j) {
        double alpha = static_cast<double>(j) / static_cast<double>(new_degree + 1);
        new_cps[cp_handler_new.Get1DIndex({j, i})] = (bezier_cps[cp_handler_old.Get1DIndex({j, i})] -
            new_cps[cp_handler_new.Get1DIndex({j - 1, i})] * alpha) * (1 / (1 - alpha));
      }
      for (int j = new_degree - 1; j >= r + 1; --j) {
        double alpha = static_cast<double>(j + 1) / static_cast<double>(new_degree + 1);
        new_cps[cp_handler_new.Get1DIndex({j, i})] = (bezier_cps[cp_handler_old.Get1DIndex({j + 1, i})] -
            new_cps[cp_handler_new.Get1DIndex({j + 1, i})] * (1 - alpha)) * (1 / alpha);
      }
      if ((new_degree + 1) % 2 != 0) {
        baf::ControlPoint P_r_L = (bezier_cps[cp_handler_old.Get1DIndex({r, i})] -
            new_cps[cp_handler_new.Get1DIndex({r - 1, i})] * alpha_r) * (1 / (1 - alpha_r));
        baf::ControlPoint P_r_R = (bezier_cps[cp_handler_old.Get1DIndex({r + 1, i})] -
            new_cps[cp_handler_new.Get1DIndex({r + 1, i})] * (1 - alpha_r_plus)) * (1 / alpha_r_plus);
        new_cps[cp_handler_new.Get1DIndex({r, i})] = (P_r_L + P_r_R) * 0.5;
        if ((P_r_L - P_r_R).GetEuclideanNorm() > max_error_bound) max_error_bound = (P_r_L - P_r_R).GetEuclideanNorm();
      } else {
        double current_max_error_bound =
            (bezier_cps[cp_handler_old.Get1DIndex({r + 1, i})] - ((new_cps[cp_handler_new.Get1DIndex({r, i})] +
                new_cps[cp_handler_new.Get1DIndex({r + 1, i})]) * 0.5)).GetEuclideanNorm();
        if (current_max_error_bound > max_error_bound) max_error_bound = current_max_error_bound;
      }
    }
    *successful = max_error_bound < tolerance;
    return new_cps;
  }

  virtual void SetNewControlPoint(baf::ControlPoint control_point, double weight,
      std::array<int, PARAMETRIC_DIMENSIONALITY> indices) = 0;

  void SetNewBezierSegmentControlPoints(const std::vector<std::vector<baf::ControlPoint>> &bezier_segments,
                                        int dimension) {
    util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> point_handler(GetPointsPerDirection());
    int width = GetDegree(dimension).Get() + 1;
    for (int i = 0; i < point_handler.GetNumberOfTotalMultiIndices(); ++i, ++point_handler) {
      int index_in_dir = point_handler[Dimension{dimension}];
      int segment = index_in_dir / (width - 1);
      int index = index_in_dir % (width - 1) + point_handler.CollapseDimension(Dimension{dimension}) * width;
      if (static_cast<size_t>(segment) > bezier_segments.size() - 1) {
        --segment;
        index = (width - 1) + point_handler.CollapseDimension(Dimension{dimension}) * width;
      }
      baf::ControlPoint cp(GetPointDim());
      double weight = bezier_segments[segment][index].GetValueForDimension(Dimension{GetPointDim()});
      for (int j = 0; j < GetPointDim(); ++j) {
        cp.SetValue(j, bezier_segments[segment][index].GetValueForDimension(Dimension{j}) / weight);
      }
      this->SetNewControlPoint(cp, weight, point_handler.GetCurrentIndex());
    }
  }

  std::shared_ptr<ParameterSpace<PARAMETRIC_DIMENSIONALITY>> parameter_space_;
};
}  // namespace splinelib::src::spl

#endif  // SRC_SPL_SPLINE_H_
