/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_B_SPLINE_H_
#define SRC_B_SPLINE_H_

#include <algorithm>
#include <functional>
#include <numeric>
#include <vector>
#include <array>

#include "control_point.h"
#include "integration_rule.h"
#include "knot_vector.h"
#include "parameter_space.h"
#include "multi_index_handler.h"

namespace spl {
template<int DIM>
class BSpline {
 public:
  BSpline(const std::array<baf::KnotVector, DIM> &knot_vector,
          std::array<int, DIM> degree,
          const std::vector<baf::ControlPoint> &control_points)
      : dim(control_points[0].GetDimension()) {
    for (int i = 0; i < DIM; ++i) {
      parameter_space_[i] = ParameterSpace(knot_vector[i], degree[i]);
    }
    for (auto &&cp : control_points) {
      for (int i = 0; i < dim; ++i) {
        control_points_.emplace_back(cp.GetValue(i));
      }
    }
  }

  std::vector<double> Evaluate(std::array<double, DIM> param_coord, const std::vector<int> &dimensions) const {
    auto basis_function_values = EvaluateAllNonZeroBasisFunctions(param_coord);
    std::vector<double> evaluated_point(dimensions.size(), 0);
    for (int i = 0; i < dimensions.size(); ++i) {
      evaluated_point[i] =
          ComputeWeightedSum(basis_function_values, ExtractControlPointValues(param_coord, dimensions[i]));
    }
    return evaluated_point;
  }

  std::vector<double> EvaluateDerivative(std::array<double, DIM> param_coord,
                                         const std::vector<int> &dimensions,
                                         std::array<int, DIM> derivative) const {
    auto basis_function_values = EvaluateAllNonZeroBasisFunctionDerivatives(param_coord, derivative);
    std::vector<double> evaluated_point(dimensions.size(), 0);
    for (int i = 0; i < dimensions.size(); ++i) {
      evaluated_point[i] =
          ComputeWeightedSum(basis_function_values, ExtractControlPointValues(param_coord, dimensions[i]));
    }
    return evaluated_point;
  }

  int GetDegree(int i) const {
    return parameter_space_[i].degree();
  }

  baf::KnotVector GetKnotVector(int i) const {
    return parameter_space_[i].knot_vector();
  }

  std::vector<elm::Element> GetElementList() const {
    return parameter_space_[0].GetElementList();
  }

  std::vector<elm::ElementIntegrationPoint> EvaluateAllElementNonZeroBasisFunctions(
      int element_number,
      const itg::IntegrationRule<1> &rule) const {
    return parameter_space_[0].EvaluateAllElementNonZeroBasisFunctions(element_number, rule);
  }

  std::vector<elm::ElementIntegrationPoint> EvaluateAllElementNonZeroBasisFunctionDerivatives(
      int element_number,
      const itg::IntegrationRule<1> &rule) const {
    return ParameterSpace2PhysicalSpace(
        parameter_space_[0].EvaluateAllElementNonZeroBasisFunctionDerivatives(element_number, rule),
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

 private:
  std::vector<double> ExtractControlPointValues(std::array<double, DIM> param_coord, int dimension) const {
    std::array<int, DIM + 1> start;
    std::array<int, DIM + 1> last;
    std::array<int, DIM + 1> total_length;
    std::array<int, DIM + 1> current;
    start[0] = dimension;
    last[0] = dimension;
    current[0] = dimension;
    total_length[0] = dim;
    int M = 1;
    for (int i = 0; i < DIM; ++i) {
      start[i + 1] = GetKnotVector(i).GetKnotSpan(param_coord[i]) - GetDegree(i);
      last[i + 1] = start[i + 1] + parameter_space_[i].degree() + 1;
      total_length[i + 1] = parameter_space_[i].knot_vector().NumberOfKnots() - parameter_space_[i].degree() - 1;
      current[i + 1] = start[i + 1];
      M *= (last[i + 1] - start[i + 1]);
    }
    util::MultiIndexHandler<DIM + 1> multiIndexHandler(total_length);
    std::vector<double> vector;
    for (int i = 0; i < M; ++i) {
      multiIndexHandler.SetIndices(current);
      vector.push_back(control_points_[multiIndexHandler.Get1DIndex()]);
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

  double ComputeWeightedSum(const std::vector<double> &basis_function_values,
                            std::vector<double> control_point_values) const {
    std::transform(basis_function_values.begin(),
                   basis_function_values.end(),
                   control_point_values.begin(),
                   control_point_values.begin(),
                   std::multiplies<double>());
    return std::accumulate(control_point_values.begin(), control_point_values.end(), 0.0, std::plus<double>());
  }

  std::vector<std::vector<double>> ParameterSpace2PhysicalSpace(std::vector<std::vector<double>> values,
                                                                int element_number,
                                                                const itg::IntegrationRule<1> &rule) const {
    elm::Element element = GetElementList()[element_number];
    for (int point = 0; point < rule.GetNumberOfIntegrationPoints(); point++) {
      std::transform(values[point].cbegin(),
                     values[point].cend(),
                     values[point].begin(),
                     std::bind2nd(std::divides<double>(),
                                  EvaluateDerivative(
                                      {ReferenceSpace2ParameterSpace(element.node(0),
                                                                     element.node(1),
                                                                     rule.coordinate(point, 0))}, {0}, {1})[0]));
    }
    return values;
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
                     std::bind2nd(std::divides<double>(),
                                  EvaluateDerivative(
                                      {ReferenceSpace2ParameterSpace(element.node(0),
                                                                     element.node(1),
                                                                     rule.coordinate(i, 0))}, {0}, {1})[0]));

      element_integration_points[i] = elm::ElementIntegrationPoint(element_non_zero_basis_functions);
    };
    return element_integration_points;
  }

  double ReferenceSpace2ParameterSpace(double upper, double lower, double point) const {
    return parameter_space_[0].ReferenceSpace2ParameterSpace(upper, lower, point);
  }

  std::array<std::vector<std::unique_ptr<baf::BasisFunction>>::const_iterator, DIM>
  CreateArrayFirstNonZeroBasisFunction(std::array<double, DIM> param_coord) const {
    std::array<std::vector<std::unique_ptr<baf::BasisFunction>>::const_iterator, DIM> first_non_zero;
    for (int i = 0; i < DIM; ++i) {
      first_non_zero[i] = parameter_space_[i].GetFirstNonZeroKnot(param_coord[i]);
    }
    return first_non_zero;
  }

  std::array<int, DIM> ArrayTotalLength() const {
    std::array<int, DIM> total_length;
    for (int i = 0; i < DIM; ++i) {
      total_length[i] = parameter_space_[i].degree() + 1;
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

  std::vector<double> EvaluateAllNonZeroBasisFunctions(std::array<double, DIM> param_coord) const {
    auto first_non_zero = this->CreateArrayFirstNonZeroBasisFunction(param_coord);
    auto total_length = this->ArrayTotalLength();
    auto M = MultiIndexHandlerShort();

    util::MultiIndexHandler<DIM> multiIndexHandler(total_length);

    std::vector<double> vector(M, 1);
    for (int i = 0; i < M; ++i) {
      for (int j = 0; j < DIM; ++j) {
        vector[i] *= (*(first_non_zero[j] + multiIndexHandler[j]))->Evaluate(param_coord[j]);
      }
      multiIndexHandler++;
    }
    return vector;
  }

  std::vector<double> EvaluateAllNonZeroBasisFunctionDerivatives(std::array<double, DIM> param_coord,
                                                                 std::array<int, DIM> derivative) const {
    auto first_non_zero = this->CreateArrayFirstNonZeroBasisFunction(param_coord);
    auto total_length = this->ArrayTotalLength();
    auto M = MultiIndexHandlerShort();

    util::MultiIndexHandler<DIM> multiIndexHandler(total_length);

    std::vector<double> vector(M, 1);
    for (int i = 0; i < M; ++i) {
      for (int j = 0; j < DIM; ++j) {
        vector[i] *= (*(first_non_zero[j] + multiIndexHandler[j]))->EvaluateDerivative(derivative[j], param_coord[j]);
      }
      multiIndexHandler++;
    }
    return vector;
  }

  std::array<ParameterSpace, DIM> parameter_space_;
  std::vector<double> control_points_;
  int dim;
};
}

#endif  // SRC_B_SPLINE_H_
