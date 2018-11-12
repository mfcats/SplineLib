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
  Spline(std::array<std::shared_ptr<baf::KnotVector>, DIM> knot_vector, std::array<Degree, DIM> degree) {
    parameter_space_ = std::make_shared<ParameterSpace<DIM>>(ParameterSpace<DIM>(knot_vector, degree));
  }
  explicit Spline(std::shared_ptr<ParameterSpace<DIM>> parameter_space) {
    parameter_space_ = parameter_space;
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

  Degree GetDegree(int i) const {
    return parameter_space_->GetDegree(i);
  }

  std::shared_ptr<baf::KnotVector> GetKnotVector(int i) const {
    return parameter_space_->GetKnotVector(i);
  }

  virtual int GetNumberOfControlPoints() = 0;
  virtual std::array<int, DIM> GetPointsPerDirection() = 0;
  virtual int GetDimension() = 0;
  virtual double GetControlPoint(std::array<int, DIM> indices, int dimension) = 0;

  virtual std::shared_ptr<spl::PhysicalSpace<DIM>> GetPhysicalSpace() const = 0;

  std::shared_ptr<spl::ParameterSpace<DIM>> GetParameterSpace() const {
    return parameter_space_;
  }

  double GetExpansion() const {
    return GetPhysicalSpace()->GetExpansion();
  }

  std::vector<double> GetControlPoints() const {
    return GetPhysicalSpace()->GetControlPoints();
  }

  std::vector<double> GetWeights() const {
    return GetPhysicalSpace()->GetWeights();
  }

  std::array<std::vector<ParamCoord>, DIM> GetKnots() const {
    return parameter_space_->GetKnots();
  }

  virtual void InsertKnot(ParamCoord knot, int dimension) {
    std::vector<double> scaling;
    KnotSpan knot_span = parameter_space_->GetKnotVector(dimension)->GetKnotSpan(knot);
    auto first = static_cast<size_t>(knot_span.get() - parameter_space_->GetDegree(dimension).get() + 1);
    auto last = static_cast<size_t>(knot_span.get());
    for (size_t i = first; i <= last; ++i) {
      ParamCoord low_knot = parameter_space_->GetKnotVector(dimension)->GetKnot(i);
      ParamCoord upper_knot =
          parameter_space_->GetKnotVector(dimension)->GetKnot(i + parameter_space_->GetDegree(dimension).get());
      scaling.emplace_back((knot.get() - low_knot.get()) / (upper_knot.get() - low_knot.get()));
    }
    this->AdjustControlPoints(scaling, static_cast<int>(first), static_cast<int>(last), dimension);
    parameter_space_->InsertKnot(knot, dimension);
  }

  virtual void AdjustControlPoints(std::vector<double> scaling, int first, int last, int dimension) = 0;

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

  std::shared_ptr<ParameterSpace<DIM>> parameter_space_;
};
}  //  namespace spl

#endif  //  SRC_SPL_SPLINE_H_
