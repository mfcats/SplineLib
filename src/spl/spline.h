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
#include "integration_rule.h"
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
  Spline(std::shared_ptr<std::array<baf::KnotVector, DIM>> knot_vector, std::array<int, DIM> degree) :
      parameter_space_(knot_vector, degree) {}
  explicit Spline(ParameterSpace<DIM> parameter_space) : parameter_space_(std::move(parameter_space)) {}

  virtual std::vector<double> Evaluate(std::array<ParamCoord, DIM> param_coord,
                                       const std::vector<int> &dimensions) const {
    ThrowIfParametricCoordinateOutsideKnotVectorRange(param_coord);

    auto first_non_zero = GetArrayOfFirstNonZeroBasisFunctions(param_coord);
    util::MultiIndexHandler<DIM> basisFunctionHandler(this->GetNumberOfBasisFunctionsToEvaluate());
    std::vector<double> evaluated_point(dimensions.size(), 0);

    for (int i = 0; i < basisFunctionHandler.Get1DLength(); ++i, basisFunctionHandler++) {
      auto indices = basisFunctionHandler.GetIndices();
      std::transform(indices.begin(), indices.end(), first_non_zero.begin(), indices.begin(), std::plus<double>());
      for (int j = 0; j < dimensions.size(); ++j) {
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
      for (int j = 0; j < dimensions.size(); ++j) {
        evaluated_point[j] += GetEvaluatedDerivativeControlPoint(param_coord, derivative, indices, dimensions[j]);
      }
    }
    return evaluated_point;
  }

  int GetDegree(int i) const {
    return parameter_space_.GetDegree(i);
  }

  std::shared_ptr<baf::KnotVector> GetKnotVector(int i) const {
    return parameter_space_[i].knot_vector();
  }

  std::vector<elm::Element> GetElementList() const {
    return parameter_space_.GetElementList(0);
  }

  std::vector<elm::ElementIntegrationPoint> EvaluateAllElementNonZeroBasisFunctions(
      int element_number,
      const itg::IntegrationRule<1> &rule) const {
    return parameter_space_.EvaluateAllElementNonZeroBasisFunctions(0, element_number, rule);
  }

  std::vector<elm::ElementIntegrationPoint> EvaluateAllElementNonZeroBasisFunctionDerivatives(
      int element_number,
      const itg::IntegrationRule<1> &rule) const {
    return ParameterSpace2PhysicalSpace(
        parameter_space_.EvaluateAllElementNonZeroBasisFunctionDerivatives(0, element_number, rule),
        element_number,
        rule);
  }

  double JacobianDeterminant(int element_number, int integration_point, const itg::IntegrationRule<1> &rule) const {
    elm::Element element = GetElementList()[element_number];
    double dx_dxi = EvaluateDerivative({ReferenceSpace2ParameterSpace(element.node(0),
                                                                      element.node(1),
                                                                      rule.coordinate(integration_point, 0))},
                                       {0},
                                       {1})[0];
    double dxi_dtildexi = (element.node(1) - element.node(0)) / 2.0;
    return dx_dxi * dxi_dtildexi;
  }

 protected:
  void ThrowIfParametricCoordinateOutsideKnotVectorRange(std::array<ParamCoord, DIM> param_coord) const {
    for (int dim = 0; dim < DIM; dim++) {
      if (!this->GetKnotVector(dim)->IsInKnotVectorRange(param_coord[dim])) {
        std::stringstream message;
        message << "The parametric coordinate " << param_coord[dim].get() << " is outside the knot vector range from "
                << GetKnotVector(dim)->GetKnot(0).get() << " to " << GetKnotVector(dim)->GetLastKnot().get() << ".";
        throw std::range_error(message.str());
      }
    }
  }

  virtual double GetEvaluatedControlPoint(std::array<ParamCoord, DIM> param_coord,
                                          std::array<int, DIM> indices,
                                          int dimension) const = 0;

  virtual double GetEvaluatedDerivativeControlPoint(std::array<ParamCoord, DIM> param_coord,
                                                    std::array<int, DIM> derivative,
                                                    std::array<int, DIM> indices,
                                                    int dimension) const = 0;

  std::vector<elm::ElementIntegrationPoint> ParameterSpace2PhysicalSpace(
      std::vector<elm::ElementIntegrationPoint> element_integration_points,
      int element_number,
      const itg::IntegrationRule<1> &rule) const {
    elm::Element element = GetElementList()[element_number];
    for (int i = 0; i < rule.GetNumberOfIntegrationPoints(); ++i) {
      std::vector<double> element_non_zero_basis_functions = element_integration_points[i].non_zero_basis_functions();
      std::transform(element_non_zero_basis_functions.cbegin(),
                     element_non_zero_basis_functions.cend(),
                     element_non_zero_basis_functions.begin(),
                     std::bind(std::divides<double>(), std::placeholders::_1,
                               EvaluateDerivative(
                                   {ReferenceSpace2ParameterSpace(element.node(0),
                                                                  element.node(1),
                                                                  rule.coordinate(i, 0))}, {0}, {1})[0]));

      element_integration_points[i] = elm::ElementIntegrationPoint(element_non_zero_basis_functions);
    }
    return element_integration_points;
  }

  ParamCoord ReferenceSpace2ParameterSpace(double upper, double lower, double point) const {
    return parameter_space_.ReferenceSpace2ParameterSpace(upper, lower, point);
  }

  std::array<int, DIM> GetArrayOfFirstNonZeroBasisFunctions(std::array<ParamCoord, DIM> param_coord) const {
    std::array<int, DIM> first_non_zero;
    for (int i = 0; i < DIM; ++i) {
      first_non_zero[i] =
          this->parameter_space_.GetKnotVector(i).GetKnotSpan(param_coord[i]) - this->parameter_space_.GetDegree(i);
    }
    return first_non_zero;
  }

  std::array<int, DIM> GetNumberOfBasisFunctionsToEvaluate() const {
    std::array<int, DIM> total_length;
    for (int i = 0; i < DIM; ++i) {
      total_length[i] = parameter_space_.GetDegree(i) + 1;
    }
    return total_length;
  }

  ParameterSpace<DIM> parameter_space_;
};
}  // namespace spl

#endif  // SRC_SPL_SPLINE_H_
