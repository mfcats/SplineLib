/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "parameter_space.h"

#include "basis_function_factory.h"
#include "element_generator.h"

spl::ParameterSpace::ParameterSpace(const baf::KnotVector &knot_vector, int degree) : degree_(degree),
                                                                                      knot_vector_(knot_vector) {
  baf::BasisFunctionFactory factory;
  basis_functions_.reserve(knot_vector_.NumberOfKnots() - degree_ - 1);
  for (uint64_t i = 0; i < (knot_vector_.NumberOfKnots() - degree_ - 1); ++i) {
    basis_functions_.emplace_back(factory.CreateDynamic(knot_vector_, i, degree_));
  }
}

std::vector<double> spl::ParameterSpace::EvaluateAllNonZeroBasisFunctions(ParamCoord param_coord) const {
  auto first_non_zero = GetFirstNonZeroKnot(param_coord);
  std::vector<double> basis_function_values(static_cast<u_int64_t >(degree_) + 1, 0.0);
  for (int i = 0; i < degree_ + 1; ++i) {
    basis_function_values[i] = (*first_non_zero)->Evaluate(param_coord);
    ++first_non_zero;
  }
  return basis_function_values;
}

std::vector<double> spl::ParameterSpace::EvaluateAllNonZeroBasisFunctionDerivatives(ParamCoord param_coord,
                                                                                    int derivative) const {
  auto first_non_zero = GetFirstNonZeroKnot(param_coord);
  std::vector<double> basis_function_values(static_cast<u_int64_t >(degree_) + 1, 0.0);
  for (int i = 0; i < degree_ + 1; ++i) {
    basis_function_values[i] = (*first_non_zero)->EvaluateDerivative(param_coord, derivative);
    ++first_non_zero;
  }
  return basis_function_values;
}

int spl::ParameterSpace::degree() const {
  return degree_;
}

baf::KnotVector spl::ParameterSpace::knot_vector() const {
  return knot_vector_;
}

std::vector<std::unique_ptr<baf::BasisFunction>>::const_iterator spl::ParameterSpace::GetFirstNonZeroKnot(
    ParamCoord param_coord) const {
  return basis_functions_.begin() + knot_vector_.GetKnotSpan(param_coord) - degree_;
}

std::vector<elm::ElementIntegrationPoint>
spl::ParameterSpace::EvaluateAllElementNonZeroBasisFunctions(int element_number,
                                                             const itg::IntegrationRule<1> &rule) const {
  elm::Element element = GetElementList()[element_number];
  std::vector<elm::ElementIntegrationPoint> element_integration_points;
  std::vector<double> non_zero_basis_functions;
  for (int i = 0; i < rule.GetNumberOfIntegrationPoints(); ++i) {
    double integration_point = ReferenceSpace2ParameterSpace(element.node(1), element.node(0), rule.coordinate(i, 0));
    non_zero_basis_functions = EvaluateAllNonZeroBasisFunctions(ParamCoord{integration_point});
    element_integration_points.emplace_back(elm::ElementIntegrationPoint(non_zero_basis_functions));
  }
  return element_integration_points;
}

std::vector<elm::ElementIntegrationPoint>
spl::ParameterSpace::EvaluateAllElementNonZeroBasisFunctionDerivatives(int element_number,
                                                                       const itg::IntegrationRule<1> &rule) const {
  elm::Element element = GetElementList()[element_number];
  std::vector<elm::ElementIntegrationPoint> element_integration_points;
  std::vector<double> non_zero_basis_function_derivatives;
  for (int i = 0; i < rule.GetNumberOfIntegrationPoints(); ++i) {
    double integration_point = ReferenceSpace2ParameterSpace(element.node(1), element.node(0), rule.coordinate(i, 0));
    non_zero_basis_function_derivatives = EvaluateAllNonZeroBasisFunctionDerivatives(ParamCoord{integration_point}, 1);
    element_integration_points.emplace_back(elm::ElementIntegrationPoint(non_zero_basis_function_derivatives));
  }
  return element_integration_points;
}

std::vector<elm::Element> spl::ParameterSpace::GetElementList() const {
  return elm::ElementGenerator(degree_, knot_vector_).GetElementList();
}

double spl::ParameterSpace::ReferenceSpace2ParameterSpace(double upper, double lower, double point) const {
  return ((upper - lower) * point + (upper + lower)) / 2.0;
}
