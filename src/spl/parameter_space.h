/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_SPL_PARAMETER_SPACE_H_
#define SRC_SPL_PARAMETER_SPACE_H_

#include <vector>

#include "basis_function.h"
#include "basis_function_factory.h"
#include "element_generator.h"
#include "element.h"
#include "integration_rule.h"
#include "knot_vector.h"
#include "element_integration_point.h"

namespace spl {
template<int DIM>
class ParameterSpace {
 public:
  ParameterSpace() = default;

  ParameterSpace(const std::array<baf::KnotVector, DIM> &knot_vector, std::array<int, DIM> degree) : degree_(degree),
                                                                                                     knot_vector_(
                                                                                                         knot_vector) {
    baf::BasisFunctionFactory factory;
    for (int i = 0; i < DIM; i++) {
      basis_functions_[i].reserve(knot_vector_[i].GetNumberOfKnots() - degree_[i] - 1);
      for (uint64_t j = 0; j < (knot_vector_[i].GetNumberOfKnots() - degree_[i] - 1); ++j) {
        basis_functions_[i].emplace_back(factory.CreateDynamic(knot_vector_[i], j, degree_[i]));
      }
    }
  }

  std::vector<double> EvaluateAllNonZeroBasisFunctions(int direction, ParamCoord param_coord) const {
    auto first_non_zero = GetFirstNonZeroKnot(direction, param_coord);
    std::vector<double> basis_function_values(static_cast<u_int64_t >(degree_[direction]) + 1, 0.0);
    for (int i = 0; i < degree_[direction] + 1; ++i) {
      basis_function_values[i] = (*first_non_zero)->Evaluate(param_coord);
      ++first_non_zero;
    }
    return basis_function_values;
  }

  std::vector<double> EvaluateAllNonZeroBasisFunctionDerivatives(int direction,
                                                                 ParamCoord param_coord,
                                                                 int derivative) const {
    auto first_non_zero = GetFirstNonZeroKnot(direction, param_coord);
    std::vector<double> basis_function_values(static_cast<u_int64_t >(degree_[direction]) + 1, 0.0);
    for (int i = 0; i < degree_[direction] + 1; ++i) {
      basis_function_values[i] = (*first_non_zero)->EvaluateDerivative(param_coord, derivative);
      ++first_non_zero;
    }
    return basis_function_values;
  }

  std::vector<std::shared_ptr<baf::BasisFunction>>::const_iterator GetFirstNonZeroKnot(int direction,
                                                                                       ParamCoord param_coord) const {
    return basis_functions_[direction].begin() + knot_vector_[direction].GetKnotSpan(param_coord) - degree_[direction];
  }

  int GetDegree(int direction) const {
    return degree_[direction];
  }

  baf::KnotVector GetKnotVector(int direction) const {
    return knot_vector_[direction];
  }

  double GetBasisFunctions(std::array<int, DIM> indices, std::array<ParamCoord, DIM> param_coord) const {
    double value = 1;
    for (int i = 0; i < DIM; ++i) {
      value *= basis_functions_[i][indices[i]]->Evaluate(param_coord[i]);
    }
    return value;
  }

  double GetBasisFunctionDerivatives(std::array<int, DIM> indices,
                                     std::array<ParamCoord, DIM> param_coord,
                                     std::array<int, DIM> derivative) const {
    double value = 1;
    for (int i = 0; i < DIM; ++i) {
      value *= basis_functions_[i][indices[i]]->EvaluateDerivative(param_coord[i], derivative[i]);
    }
    return value;
  }

  std::vector<elm::ElementIntegrationPoint>
  EvaluateAllElementNonZeroBasisFunctions(int direction,
                                          int element_number,
                                          const itg::IntegrationRule<1> &rule) const {
    elm::Element element = GetElementList(direction)[element_number];
    std::vector<elm::ElementIntegrationPoint> element_integration_points;
    std::vector<double> non_zero_basis_functions;
    for (int i = 0; i < rule.GetNumberOfIntegrationPoints(); ++i) {
      ParamCoord
          integration_point = ReferenceSpace2ParameterSpace(element.node(1), element.node(0), rule.coordinate(i, 0));
      non_zero_basis_functions = EvaluateAllNonZeroBasisFunctions(direction, ParamCoord{integration_point});
      element_integration_points.emplace_back(elm::ElementIntegrationPoint(non_zero_basis_functions));
    }
    return element_integration_points;
  }

  std::vector<elm::ElementIntegrationPoint>
  EvaluateAllElementNonZeroBasisFunctionDerivatives(int direction,
                                                    int element_number,
                                                    const itg::IntegrationRule<1> &rule) const {
    elm::Element element = GetElementList(direction)[element_number];
    std::vector<elm::ElementIntegrationPoint> element_integration_points;
    std::vector<double> non_zero_basis_function_derivatives;
    for (int i = 0; i < rule.GetNumberOfIntegrationPoints(); ++i) {
      ParamCoord
          integration_point = ReferenceSpace2ParameterSpace(element.node(1), element.node(0), rule.coordinate(i, 0));
      non_zero_basis_function_derivatives =
          EvaluateAllNonZeroBasisFunctionDerivatives(direction, ParamCoord{integration_point}, 1);
      element_integration_points.emplace_back(elm::ElementIntegrationPoint(non_zero_basis_function_derivatives));
    }
    return element_integration_points;
  }

  std::vector<elm::Element> GetElementList(int direction) const {
    return elm::ElementGenerator(degree_[direction], knot_vector_[direction]).GetElementList();
  }

  ParamCoord ReferenceSpace2ParameterSpace(double upper, double lower, double point) const {
    return ParamCoord{((upper - lower) * point + (upper + lower)) / 2.0};
  }

 private:
  std::array<int, DIM> degree_;
  std::array<baf::KnotVector, DIM> knot_vector_;
  std::array<std::vector<std::shared_ptr<baf::BasisFunction>>, DIM> basis_functions_;
};
}  // namespace spl

#endif  // SRC_SPL_PARAMETER_SPACE_H_
