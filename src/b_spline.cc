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
#include <numeric>

#include "basis_function_factory.h"

BSpline::BSpline(const KnotVector &knot_vector, Degree degree,
                 const std::vector<ControlPoint> &control_points) :
    parameter_space_(ParameterSpace(knot_vector, degree)) {
  for (int control_point = 0; control_point < control_points.size(); control_point++) {
    control_points_.emplace_back(std::vector<double>(0));
    for (int dimension = 0; dimension < control_points[0].dimension(); dimension++) {
      control_points_[control_point].emplace_back(control_points[control_point].GetValue(dimension));
    }
  }
}

std::vector<double> BSpline::Evaluate(double param_coord, std::vector<Dimension> dimensions) const {
  auto basis_function_values = parameter_space_.EvaluateAllNonZeroBasisFunctions(param_coord);
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
  auto basis_function_values =
      parameter_space_.EvaluateAllNonZeroBasisFunctionDerivatives(param_coord, derivative);
  std::vector<double> evaluated_point(dimensions.size(), 0);
  for (int i = 0; i < dimensions.size(); ++i) {
    evaluated_point[i] = ComputeWeightedSum(basis_function_values,
                                            ExtractControlPointValues(param_coord, dimensions[i]));
  }
  return evaluated_point;
}

int BSpline::degree() const {
  return parameter_space_.degree();
}

KnotVector BSpline::knot_vector() const {
  return parameter_space_.knot_vector();
}

std::vector<double> BSpline::ExtractControlPointValues(double param_coord, int dimension) const {
  std::vector<double> control_point_values(static_cast<uint64_t>(degree() + 1), 0.0);
  auto control_points = control_points_.begin() + knot_vector().GetKnotSpan(param_coord) - degree();
  for (int i = 0; i < degree() + 1; ++i) {
    control_point_values[i] = (*control_points)[dimension];
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
