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

  Spline(std::shared_ptr<std::array<baf::KnotVector, DIM>> knot_vector,
         std::array<int, DIM> degree,
         const std::vector<baf::ControlPoint> &control_points) : parameter_space_(*knot_vector, degree) {
    std::array<int, DIM> number_of_points;
    for (int i = 0; i < DIM; ++i) {
      number_of_points[i] = (*knot_vector)[i].GetNumberOfKnots() - degree[i] - 1;
    }
    physical_space_ = PhysicalSpace<DIM>(control_points, number_of_points);
  }

  Spline(ParameterSpace<DIM> parameter_space, PhysicalSpace<DIM> physical_space) : parameter_space_(
      std::move(parameter_space)), physical_space_(physical_space) {}

  std::vector<double> Evaluate(std::array<ParamCoord, DIM> param_coord, const std::vector<int> &dimensions) const {
    ThrowIfParametricCoordinateOutsideKnotVectorRange(param_coord);
    auto basis_function_values = EvaluateAllNonZeroBasisFunctions(param_coord);
    std::vector<double> evaluated_point(dimensions.size(), 0);
    for (int i = 0; i < dimensions.size(); ++i) {
      evaluated_point[i] =
          ComputeWeightedSum(basis_function_values, ExtractControlPointValues(param_coord, dimensions[i]));
    }
    return evaluated_point;
  }

  virtual std::vector<double> EvaluateDerivative(std::array<ParamCoord, DIM> param_coord,
                                                 const std::vector<int> &dimensions,
                                                 std::array<int, DIM> derivative) const = 0;

  int GetDegree(int i) const {
    return parameter_space_.degree(i);
  }

  baf::KnotVector GetKnotVector(int i) const {
    return parameter_space_.knot_vector(i);
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

  std::vector<double> ExtractControlPointValues(std::array<ParamCoord, DIM> param_coord, int dimension) const {
    std::array<int, DIM + 1> start;
    std::array<int, DIM + 1> last;
    std::array<int, DIM + 1> total_length;
    std::array<int, DIM + 1> current;
    start[0] = dimension;
    last[0] = dimension;
    current[0] = dimension;
    total_length[0] = physical_space_.GetDimension();
    int M = 1;
    for (int i = 0; i < DIM; ++i) {
      start[i + 1] = GetKnotVector(i).GetKnotSpan(param_coord[i]) - GetDegree(i);
      last[i + 1] = start[i + 1] + parameter_space_.degree(i) + 1;
      total_length[i + 1] = parameter_space_.knot_vector(i).GetNumberOfKnots() - parameter_space_.degree(i) - 1;
      current[i + 1] = start[i + 1];
      M *= (last[i + 1] - start[i + 1]);
    }
    util::MultiIndexHandler<DIM + 1> multiIndexHandler(total_length);
    std::vector<double> vector;
    for (int i = 0; i < M; ++i) {
      multiIndexHandler.SetIndices(current);
      vector.push_back(physical_space_.GetPoint(multiIndexHandler.Get1DIndex()));
      for (int i = 0; i < DIM; ++i) {
        if (current[i + 1] == last[i + 1] - 1) {
          current[i + 1] = start[i + 1];
        } else {
          current[i + 1]++;
          break;
        }
      }
    }
    return vector;
  }

 protected:
  void ThrowIfParametricCoordinateOutsideKnotVectorRange(std::array<ParamCoord, DIM> param_coord) const {
    for (int dim = 0; dim < DIM; dim++) {
      if (!this->GetKnotVector(dim).IsInKnotVectorRange(param_coord[dim])) {
        std::stringstream message;
        message << "The parametric coordinate " << param_coord[dim].get() << " is outside the knot vector range from "
                << GetKnotVector(dim).GetKnot(0).get() << " to " << GetKnotVector(dim).GetLastKnot().get() << ".";
        throw std::range_error(message.str());
      }
    }
  }

  double ComputeWeightedSum(const std::vector<double> &basis_function_values,
                            std::vector<double> control_point_values) const {
    std::transform(basis_function_values.begin(),
                   basis_function_values.end(),
                   control_point_values.begin(),
                   control_point_values.begin(),
                   std::multiplies<double>());
    return std::accumulate(control_point_values.begin(), control_point_values.end(), 0.0, std::plus<double>());
  }

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

  std::array<std::vector<std::shared_ptr<baf::BasisFunction>>::const_iterator, DIM>
  CreateArrayFirstNonZeroBasisFunction(std::array<ParamCoord, DIM> param_coord) const {
    std::array<std::vector<std::shared_ptr<baf::BasisFunction>>::const_iterator, DIM> first_non_zero;
    for (int i = 0; i < DIM; ++i) {
      first_non_zero[i] = parameter_space_.GetFirstNonZeroKnot(i, param_coord[i]);
    }
    return first_non_zero;
  }

  std::array<int, DIM> ArrayTotalLength() const {
    std::array<int, DIM> total_length;
    for (int i = 0; i < DIM; ++i) {
      total_length[i] = parameter_space_.degree(i) + 1;
    }
    return total_length;
  }

  int MultiIndexHandlerShort() const {
    int M = 1;
    std::array<int, DIM> total_length = ArrayTotalLength();
    for (int i = 0; i < DIM; ++i) {
      M *= total_length[i];
    }
    return M;
  }

  virtual std::vector<double> EvaluateAllNonZeroBasisFunctions(std::array<ParamCoord, DIM> param_coord) const = 0;

  ParameterSpace<DIM> parameter_space_;
  PhysicalSpace<DIM> physical_space_;
};
}  // namespace spl

#endif  // SRC_SPL_SPLINE_H_
