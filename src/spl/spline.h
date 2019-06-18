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
#include <iostream>
#include <algorithm>
#include <array>
#include <functional>
#include <numeric>
#include <sstream>
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
  explicit Spline(std::shared_ptr<ParameterSpace < DIM>>
parameter_space) {
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
      std::transform(indices.begin(), indices.end(), first_non_zero.begin(), indices.begin(), std::plus<double>());
      for (auto j = 0u; j < dimensions.size(); ++j) {
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
      std::transform(indices.begin(), indices.end(), first_non_zero.begin(), indices.begin(), std::plus<double>());
      for (auto j = 0u; j < dimensions.size(); ++j) {
        evaluated_point[j] += GetEvaluatedDerivativeControlPoint(param_coord, derivative, indices, dimensions[j]);
      }
    }
    return evaluated_point;
  }

  std::vector<double> EvaluateAllNonZeroBasisFunctions(int direction, ParamCoord param_coord) const {
    return parameter_space_->EvaluateAllNonZeroBasisFunctions(direction, param_coord);
  }

  std::vector<double> EvaluateAllNonZeroBasisFunctionDerivatives(int direction,
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

  std::array<int, DIM> GetPointsPerDirection() const {
    return GetPhysicalSpace()->GetPointsPerDirection();
  }

  int GetPointDim() const {
    return GetPhysicalSpace()->GetDimension();
  }

  double GetControlPoint(std::array<int, DIM> indices, int dimension) const {
    return GetPhysicalSpace()->GetControlPoint(indices).GetValue(dimension);
  }

  baf::ControlPoint GetControlPoint(std::array<int, DIM> indices) const {
    return GetPhysicalSpace()->GetControlPoint(indices);
  }

  double GetExpansion() const {
    return GetPhysicalSpace()->GetExpansion();
  }

  double GetWeight(std::array<int, DIM> indices) const {
    return GetPhysicalSpace()->GetWeight(indices);
  }

  void InsertKnot(ParamCoord knot, int dimension, size_t multiplicity = 1) {
    KnotSpan knot_span = GetKnotVector(dimension)->GetKnotSpan(knot);
    Degree degree = GetDegree(dimension);
    for (size_t i = 1; i <= multiplicity; ++i) {
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

  void RefineKnots(std::vector<ParamCoord> new_knots, int dimension) {
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

  void ElevateDegree(int dimension) {
    std::vector<double> alpha = ComputeBezierDegreeElevationCoeffients(dimension);
    auto diff = ProduceBezierSegments(dimension);

    std::vector<std::vector<baf::ControlPoint>>
        all_new_bezier_cps(GetKnotVector(dimension)->GetNumberOfDifferentKnots() - 1);
    for (int i = static_cast<int>(GetKnotVector(dimension)->GetNumberOfDifferentKnots() - 2); i >= 0; --i) {
      std::vector<double> bezier_cps = GetBezierSegment(dimension, i);
      all_new_bezier_cps[i] = DegreeElevateBezierSegment(bezier_cps, alpha, dimension);
    }
    InsertKnot(GetKnotVector(dimension)->GetLastKnot(), dimension);
    parameter_space_->InsertKnot(GetKnotVector(dimension)->GetKnot(0), dimension);
    for (int i = static_cast<int>(GetKnotVector(dimension)->GetNumberOfDifferentKnots() - 2); i >= 0; --i) {
      SetNewBezierSegmentControlPoints(all_new_bezier_cps[i], dimension, i);
    }
    parameter_space_->ElevateDegree(dimension);
    RemoveBezierKnots(diff, dimension);
  }

  void ReduceDegree(int dimension) {
    // Insert knots in order to produce bezier segments.
    auto diff = ProduceBezierSegments(dimension);
    // Vector to hold the bezier segment control points for all bezier segments.
    std::vector<std::vector<baf::ControlPoint>>
        all_new_bezier_cps(GetKnotVector(dimension)->GetNumberOfDifferentKnots() - 1);
    // Loop backwards through the created bezier segments.
    for (int i = static_cast<int>(GetKnotVector(dimension)->GetNumberOfDifferentKnots() - 2); i >= 0; --i) {
      // Store the control point coordinates of the current bezier segment (x1, y1, z1, x2, y2, z2, ...).
      std::vector<double> bezier_cps = GetBezierSegment(dimension, i);
      all_new_bezier_cps[i] = DegreeReduceBezierSegment(bezier_cps, dimension);
    }
    // Assemble the computed control points into the physical space of the spline.
    for (int i = static_cast<int>(GetKnotVector(dimension)->GetNumberOfDifferentKnots() - 2); i >= 0; --i) {
      // TODO: Create one SetNewBezierSegmentControlPoints function for degree elevation and reduction that also adds
      // TODO: or removes the necessary control points.
      SetNewBezierSegmentControlPoints2(all_new_bezier_cps[i], dimension, i);
    }
    GetPhysicalSpace()->RemoveControlPoints(GetKnotVector(dimension)->GetNumberOfDifferentKnots() - 1);
    // Remove one entry of the first and last knot.
    parameter_space_->RemoveKnot(GetKnotVector(dimension)->GetKnot(0), dimension);
    parameter_space_->RemoveKnot(GetKnotVector(dimension)->GetLastKnot(), dimension);
    // Reduce degree of the dimension by one and create new basis functions.
    parameter_space_->ReduceDegree(dimension);
    // Remove the knots to get rid of the decomposition into bezier segments (remove one knot more than inserted).
    std::transform(diff.begin(), diff.end(), diff.begin(), std::bind2nd(std::plus<int>(), 2));
    // TODO: Segfault occurs when executing the function below.
    RemoveBezierKnots(diff, dimension);
  }

  std::vector<baf::ControlPoint> DegreeReduceBezierSegment(std::vector<double> bezier_cps, int dimension) const {
    uint64_t r = (GetDegree(dimension).get() - 1) / 2;
    // New control points.
    std::vector<std::vector<double>> cps(static_cast<unsigned long>(GetDegree(dimension).get()));
    // Vector of length sdim (number of space dimensions).
    std::vector<double> coord(static_cast<unsigned long>(GetPointDim()));
    // First control point of the degree reduced bezier segment is equal to the original first control point.
    for (uint64_t j = 0; j < GetPointDim(); ++j) {
      coord[j] = bezier_cps[j];
    }
    cps[0] = coord;
    // Last control point of the degree reduced bezier segment is equal to the original last control point.
    for (uint64_t j = 0; j < GetPointDim(); ++j) {
      coord[j] = bezier_cps[GetDegree(dimension).get() * GetPointDim() + j];
    }
    cps[GetDegree(dimension).get() - 1] = coord;
    // Loop through the new control points P_i for i = 1,...,r-1.
    for (uint64_t i = 1; i <= r - 1; ++i) {
      // Compute the degree elevation coeffiecient alpha_i for the current control point.
      double alpha_i = double(i) / double(GetDegree(dimension).get());
      // Loop through the physical dimensions of the current new control point.
      for (uint64_t j = 0; j < GetPointDim(); ++j) {
        // Compute the new control point P_i as a function of P_i-1, the old control point Q_i and the degree elevation
        // coefficient alpha_i.
        coord[j] = (bezier_cps[i * GetPointDim() + j] - alpha_i * cps[i - 1][j]) / (1 - alpha_i);
      }
      cps[i] = coord;
    }
    // Compute the new control point P_r depending on whether the degree p is even or odd.
    if (GetDegree(dimension).get() % 2 == 0) {
      double alpha_i = double(r) / double(GetDegree(dimension).get());
      for (uint64_t j = 0; j < GetPointDim(); ++j) {
        coord[j] = (bezier_cps[r * GetPointDim() + j] - alpha_i * cps[r - 1][j]) / (1 - alpha_i);
      }
    } else {
      double alpha_r_L = double(r) / double(GetDegree(dimension).get());
      double alpha_r_R = double(r + 1) / double(GetDegree(dimension).get());
      for (uint64_t j = 0; j < GetPointDim(); ++j) {
        double P_r_L = (bezier_cps[r * GetPointDim() + j] - alpha_r_L * cps[r - 1][j]) / (1 - alpha_r_L);
        double P_r_R = (bezier_cps[(r + 1) * GetPointDim() + j] - (1 - alpha_r_R) * cps[(r + 1) - 1][j])
                       / (1 - alpha_r_R);
        coord[j] = 0.5 * (P_r_L + P_r_R);
      }
    }
    cps[r] = coord;
    // Loop though the remaining new control points P_i for i = p-2,...,r+1.
    for (uint64_t i = r + 1; i <= GetDegree(dimension).get() - 2; ++i) {
      // Compute the degree elevation coeffiecient alpha_i for the current control point.
      double alpha_i = double(i + 1) / double(GetDegree(dimension).get());

      // std::cout << std::endl << "alpha_i: " << alpha_i; = 0.75

      // Loop through the physical dimensions of the current new control point.
      for (uint64_t j = 0; j < GetPointDim(); ++j) {
        // Compute the new control point P_i as a function of P_i-1, the old control point Q_i and the degree elevation
        // coefficient alpha_i.
        coord[j] = (bezier_cps[(i + 1) * GetPointDim() + j] - (1 - alpha_i) * cps[i + 1][j]) / alpha_i;
      }
      cps[i] = coord;
    }
    // Create vector of baf::ControlPoint out of the cps vector.
    std::vector<baf::ControlPoint> cps_;
    for (uint64_t i = 0; i < cps.size(); ++i) {
      cps_.emplace_back(baf::ControlPoint(cps[i]));
    }
    return cps_;
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

  std::vector<double> ComputeBezierDegreeElevationCoeffients(int dimension) const {
    std::vector<double> alpha;
    for (int i = 0; i < GetDegree(dimension).get() + 2; ++i) {
      alpha.push_back(static_cast<double>(i) / (GetDegree(dimension).get() + 1));
    }
    return alpha;
  }

  std::vector<int> ProduceBezierSegments(int dimension) {
    size_t position = GetKnotVector(dimension)->GetNumberOfKnots() - GetDegree(dimension).get() - 2;
    std::vector<int> diff(GetKnotVector(dimension)->GetNumberOfDifferentKnots() - 2, GetDegree(dimension).get());
    for (size_t i = GetKnotVector(dimension)->GetNumberOfDifferentKnots() - 2; i > 0; --i, --position) {
      while (GetKnotVector(dimension)->GetKnot(position) == GetKnotVector(dimension)->GetKnot(position - 1)) {
        --position;
        --diff[i - 1];
      }
      ParamCoord cur_knot = GetKnotVector(dimension)->GetKnot(position);
      InsertKnot(cur_knot, dimension, static_cast<size_t>(diff[i - 1]));
    }
    return diff;
  }

  void RemoveBezierKnots(std::vector<int> diff, int dimension) {
    size_t position = GetKnotVector(dimension)->GetNumberOfKnots() - GetDegree(dimension).get() - 2;
    for (size_t i = GetKnotVector(dimension)->GetNumberOfDifferentKnots() - 2; i > 0; --i, --position) {
      while (GetKnotVector(dimension)->GetKnot(position) == GetKnotVector(dimension)->GetKnot(position - 1)) {
        --position;
      }
      ParamCoord cur_knot = GetKnotVector(dimension)->GetKnot(position);
      RemoveKnot(cur_knot, dimension, util::NumericSettings<double>::kEpsilon(), static_cast<size_t>(diff[i - 1] - 1));
    }
  }

  std::vector<double> GetBezierSegment(int dimension, int segment) const {
    int first = segment * (GetDegree(dimension).get() + 1);
    std::vector<double> bezier_cps(static_cast<unsigned long>(GetPointDim() * (GetDegree(dimension).get() + 1)));
    for (int i = 0; i < GetDegree(dimension).get() + 1; ++i) {
      for (int j = 0; j < GetPointDim(); ++j) {
        bezier_cps[i * GetPointDim() + j] = GetControlPoint({first + i}, j);
      }
    }
    return bezier_cps;
  }

  std::vector<baf::ControlPoint> DegreeElevateBezierSegment(std::vector<double> bezier_cps,
                                                            std::vector<double> alpha,
                                                            int dimension) const {
    std::vector<baf::ControlPoint> cps;
    std::vector<double> coord(static_cast<unsigned long>(GetPointDim()));
    for (int i = 0; i < GetDegree(dimension).get() + 1; ++i) {
      for (int j = 0; j < GetPointDim(); ++j) {
        coord[j] =
            (1 - alpha[i]) * bezier_cps[i * GetPointDim() + j] + alpha[i] * bezier_cps[(i - 1) * GetPointDim() + j];
      }
      cps.emplace_back(baf::ControlPoint(coord));
    }
    for (int j = 0; j < GetPointDim(); ++j) {
      coord[j] = bezier_cps[GetDegree(dimension).get() * GetPointDim() + j];
    }
    cps.emplace_back(baf::ControlPoint(coord));
    return cps;
  }

  void SetNewBezierSegmentControlPoints(std::vector<baf::ControlPoint> new_bezier_cps, int dimension, int segment) {
    int first = segment * (GetDegree(dimension).get() + 1);
    for (int i = 0; i < GetDegree(dimension).get() + 1; ++i) {
      GetPhysicalSpace()->SetControlPoint({i + first}, new_bezier_cps[i]);
    }
  }

  void SetNewBezierSegmentControlPoints2(std::vector<baf::ControlPoint> new_bezier_cps, int dimension, int segment) {
    int first = segment * (GetDegree(dimension).get());
    for (int i = 0; i < GetDegree(dimension).get(); ++i) {
      GetPhysicalSpace()->SetControlPoint({i + first}, new_bezier_cps[i]);
    }
  }

  std::shared_ptr<ParameterSpace < DIM>> parameter_space_;
};
}  //  namespace spl

#endif  //  SRC_SPL_SPLINE_H_
