/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#include "b_spline.h"

#include <algorithm>
#include <functional>
#include <numeric>

BSpline::BSpline(const KnotVector &knot_vector, int degree, const std::vector<ControlPoint> &control_points)
    : parameter_space_(ParameterSpace(knot_vector, degree)), spline_dimension_(control_points[0].GetDimension()) {
  for (auto &&cp : control_points) {
    for (int i = 0; i < spline_dimension_; ++i) {
      control_points_.emplace_back(cp.GetValue(i));
    }
  }
}

std::vector<double> BSpline::Evaluate(double param_coord, const std::vector<int> &dimensions) const {
  auto basis_function_values = parameter_space_.EvaluateAllNonZeroBasisFunctions(param_coord);
  std::vector<double> evaluated_point(dimensions.size(), 0);
  for (int i = 0; i < dimensions.size(); ++i) {
    evaluated_point[i] =
        ComputeWeightedSum(basis_function_values, ExtractControlPointValues(param_coord, dimensions[i]));
  }
  return evaluated_point;
}

std::vector<double> BSpline::EvaluateDerivative(double param_coord,
                                                const std::vector<int> &dimensions,
                                                int derivative) const {
  auto basis_function_values = parameter_space_.EvaluateAllNonZeroBasisFunctionDerivatives(param_coord, derivative);
  std::vector<double> evaluated_point(dimensions.size(), 0);
  for (int i = 0; i < dimensions.size(); ++i) {
    evaluated_point[i] =
        ComputeWeightedSum(basis_function_values, ExtractControlPointValues(param_coord, dimensions[i]));
  }
  return evaluated_point;
}

std::vector<Element> BSpline::GetElementList() const {
  return parameter_space_.GetElementList();
}

std::vector<std::vector<double>> BSpline::EvaluateAllElementNonZeroBasisFunctions(
    int element_number,
    const IntegrationRule<1> &rule) const {
  return parameter_space_.EvaluateAllElementNonZeroBasisFunctions(element_number, rule);
}

std::vector<std::vector<double>> BSpline::EvaluateAllElementNonZeroBasisFunctionDerivatives(
    int element_number,
    const IntegrationRule<1> &rule) const {
  return TransformToPhysicalSpace(
      parameter_space_.EvaluateAllElementNonZeroBasisFunctionDerivatives(element_number, rule),
      element_number,
      rule);
}

double BSpline::JacobianDeterminant(int element_number, int integration_point, const IntegrationRule<1> &rule) const {
  Element element = GetElementList()[element_number];
  double dx_dxi = EvaluateDerivative(TransformToParameterSpace(element.node(0),
                                                               element.node(1),
                                                               rule.coordinate(integration_point, 0)), {0}, 1)[0];
  double dxi_dtildexi = (element.node(1) - element.node(0))/2.0;
  return dx_dxi*dxi_dtildexi;
}

int BSpline::GetDegree() const {
  return parameter_space_.degree();
}

KnotVector BSpline::GetKnotVector() const {
  return parameter_space_.knot_vector();
}

std::vector<double> BSpline::ExtractControlPointValues(double param_coord, int dimension) const {
  std::vector<double> control_point_values(static_cast<uint64_t>(GetDegree() + 1), 0.0);
  auto control_points =
      control_points_.begin() + (GetKnotVector().GetKnotSpan(param_coord) - GetDegree())*spline_dimension_ + dimension;
  for (int i = 0; i < GetDegree() + 1; ++i) {
    control_point_values[i] = *control_points;
    control_points += spline_dimension_;
  }
  return control_point_values;
}

double BSpline::ComputeWeightedSum(const std::vector<double> &basis_function_values,
                                   std::vector<double> control_point_values) const {
  std::transform(basis_function_values.begin(),
                 basis_function_values.end(),
                 control_point_values.begin(),
                 control_point_values.begin(),
                 std::multiplies<double>());
  return std::accumulate(control_point_values.begin(), control_point_values.end(), 0.0, std::plus<double>());
}

std::vector<std::vector<double>> BSpline::TransformToPhysicalSpace(std::vector<std::vector<double>> values,
                                                                   int element_number,
                                                                   const IntegrationRule<1> &rule) const {
  Element element = GetElementList()[element_number];
  for (int point = 0; point < rule.GetNumberOfIntegrationPoints(); point++) {
    std::transform(values[point].cbegin(),
                   values[point].cend(),
                   values[point].begin(),
                   std::bind2nd(std::divides<double>(),
                                EvaluateDerivative(TransformToParameterSpace(element.node(0),
                                                                             element.node(1),
                                                                             rule.coordinate(point, 0)), {0}, 1)[0]));
  }
  return values;
}

double BSpline::TransformToParameterSpace(double upper, double lower, double point) const {
  return parameter_space_.TransformToParameterSpace(upper, lower, point);
}
