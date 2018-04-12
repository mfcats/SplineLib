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

#include "basis_function_factory.h"

BSpline::BSpline(const KnotVector &knot_vector,
                 Degree degree,
                 const std::vector<ControlPoint> &control_points) :
    degree_(degree), control_points_(control_points), knot_vector_(knot_vector) {
  BasisFunctionFactory factory;
  basis_functions_.reserve(control_points_.size());
  for (uint64_t i = 0; i < control_points_.size(); ++i) {
    basis_functions_.emplace_back(factory.CreateDynamic(knot_vector_, i, degree_));
  }
}

std::vector<double> BSpline::Evaluate(double param_coord, std::vector<Dimension> dimensions) const {
  auto basis_function_values = EvaluateAllNonZeroBasisFunctions(param_coord);
  std::vector<double> evaluated_point(dimensions.size(), 0);
  for (int i = 0; i < dimensions.size(); ++i) {
    evaluated_point[i] = ComputeWeightedSum(basis_function_values,
                                            ExtractControlPointValues(param_coord, dimensions[i]));
  }
  return evaluated_point;
}

std::vector<double> BSpline::EvaluateDerivative(double param_coord,
                                                std::vector<int> dimensions,
                                                int derivative) const {
  auto basis_function_values = EvaluateAllNonZeroBasisFunctionDerivatives(param_coord, derivative);
  std::vector<double> evaluated_point(dimensions.size(), 0);
  for (int i = 0; i < dimensions.size(); ++i) {
    evaluated_point[i] = ComputeWeightedSum(basis_function_values,
                                            ExtractControlPointValues(param_coord, dimensions[i]));
  }
  return evaluated_point;
}

std::vector<double> BSpline::EvaluateAllNonZeroBasisFunctions(double param_coord) const {
  auto first_non_zero = basis_functions_.begin() + knot_vector_.GetKnotSpan(param_coord) - degree_;
  std::vector<double> basis_function_values(degree_ + 1, 0);
  for (int i = 0; i < degree_ + 1; ++i) {
    basis_function_values[i] = (*first_non_zero)->Evaluate(param_coord);
    ++first_non_zero;
  }
  return basis_function_values;
}

std::vector<double> BSpline::EvaluateAllNonZeroBasisFunctionDerivatives(double param_coord,
                                                                        int derivative) const {
  auto first_non_zero = basis_functions_.begin() + knot_vector_.GetKnotSpan(param_coord) - degree_;
  std::vector<double> basis_function_values(degree_ + 1, 0);
  for (int i = 0; i < degree_ + 1; ++i) {
    basis_function_values[i] = (*first_non_zero)->EvaluateDerivative(derivative, param_coord);
    ++first_non_zero;
  }
  return basis_function_values;
}

std::vector<double> BSpline::ExtractControlPointValues(double param_coord, int dimension) const {
  std::vector<double> control_point_values(static_cast<uint64_t>(degree_ + 1), 0.0);
  auto control_points = control_points_.begin() + knot_vector_.GetKnotSpan(param_coord) - degree_;
  for (int i = 0; i < degree_ + 1; ++i) {
    control_point_values[i] = control_points->GetValue(dimension);
    ++control_points;
  }
  return control_point_values;
}

double BSpline::ComputeWeightedSum(std::vector<double> basis_function_values,
                                   std::vector<double> control_point_values) const {
  std::transform(basis_function_values.begin(),
                 basis_function_values.end(),
                 control_point_values.begin(),
                 control_point_values.begin(),
                 std::multiplies<double>());
  return std::accumulate(control_point_values.begin(),
                         control_point_values.end(),
                         0.0,
                         std::plus<double>());
}
