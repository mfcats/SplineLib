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
    // compute Bezier coefficients
    std::vector<double> alpha;
    for (double i = 0; i < GetDegree(dimension).get() + 2; ++i) {
      alpha.push_back(i / (GetDegree(dimension).get() + 1));
    }
    // b: insert knot 0.3 twice
    InsertKnot(ParamCoord{0.3}, 0, 2);
    // c: degree elevate Bezier segment [0, 0.3]
    std::array<double, 4 * 2> bezier_cp;
    std::array<double, 5 * 2> new_bezier_cp;
    std::vector<baf::ControlPoint> cps;
    for (int i = 0; i < GetDegree(dimension).get() + 1; ++i) {
      for (int j = 0; j < GetPointDim(); ++j) {
        bezier_cp[i * GetPointDim() + j] = GetControlPoint({i}, j);
      }
    }
    std::vector<double> coord(GetPointDim());
    for (int i = 0; i < GetDegree(dimension).get() + 1; ++i) {
      for (int j = 0; j < GetPointDim(); ++j) {
        double
            n = (1 - alpha[i]) * bezier_cp[i * GetPointDim() + j] + alpha[i] * bezier_cp[(i - 1) * GetPointDim() + j];
        new_bezier_cp[i * GetPointDim() + j] = n;
        coord[j] = n;
      }
      cps.emplace_back(baf::ControlPoint(coord));
    }
    for (int j = 0; j < GetPointDim(); ++j) {
      new_bezier_cp[(GetDegree(dimension).get() + 1) * GetPointDim() + j] =
          bezier_cp[GetDegree(dimension).get() * GetPointDim() + j];
      coord[j] = new_bezier_cp[(GetDegree(dimension).get() + 1) * GetPointDim() + j];
    }
    cps.emplace_back(baf::ControlPoint(coord));
    GetPhysicalSpace()->AddControlPoints(1);
    for (int i = GetNumberOfControlPoints(); i > 0; --i) {
      GetPhysicalSpace()->SetControlPoint({i}, GetPhysicalSpace()->GetControlPoint({i - 1}));
    }
    for (int i = 0; i < GetDegree(dimension).get() + 2; ++i) {
      GetPhysicalSpace()->SetControlPoint({i}, cps[i]);
    }
    GetPhysicalSpace()->IncrementNumberOfPoints(dimension);
    // d: insert knot 0.6 twice
    InsertKnot(ParamCoord{0.6}, 0, 2);
    int offset = 4;
    // e: degree elevate Bezier segment [0.3, 0.6]
    cps.erase(cps.begin(), cps.end());
    std::vector<baf::ControlPoint> cps2;
    for (int i = 0; i < GetDegree(dimension).get() + 1; ++i) {
      for (int j = 0; j < GetPointDim(); ++j) {
        bezier_cp[i * GetPointDim() + j] = GetControlPoint({i + offset}, j);
      }
    }
    for (int i = 0; i < GetDegree(dimension).get() + 1; ++i) {
      for (int j = 0; j < GetPointDim(); ++j) {
        double
            n = (1 - alpha[i]) * bezier_cp[i * GetPointDim() + j] + alpha[i] * bezier_cp[(i - 1) * GetPointDim() + j];
        new_bezier_cp[i * GetPointDim() + j] = n;
        coord[j] = n;
      }
      cps2.emplace_back(baf::ControlPoint(coord));
    }
    for (int j = 0; j < GetPointDim(); ++j) {
      new_bezier_cp[(GetDegree(dimension).get() + 1) * GetPointDim() + j] =
          bezier_cp[GetDegree(dimension).get() * GetPointDim() + j];
      coord[j] = new_bezier_cp[(GetDegree(dimension).get() + 1) * GetPointDim() + j];
    }
    cps2.emplace_back(baf::ControlPoint(coord));
    GetPhysicalSpace()->AddControlPoints(1);
    GetPhysicalSpace()->IncrementNumberOfPoints(dimension);
    for (int i = GetNumberOfControlPoints(); i > offset; --i) {
      GetPhysicalSpace()->SetControlPoint({i}, GetPhysicalSpace()->GetControlPoint({i - 1}));
    }
    for (int i = 0; i < GetDegree(dimension).get() + 2; ++i) {
      GetPhysicalSpace()->SetControlPoint({i + offset}, cps2[i]);
    }
    GetPhysicalSpace()->IncrementNumberOfPoints(dimension);
    // b: remove knot 0.3 twice
//    RemoveKnot(ParamCoord{0.3}, 0, 0.1, 2);
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

  std::shared_ptr<ParameterSpace < DIM>> parameter_space_;
};
}  //  namespace spl

#endif  //  SRC_SPL_SPLINE_H_
