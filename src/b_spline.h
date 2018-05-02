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
#include "parameter_space.h"

template <int DIM>
class BSpline {
 public:
  BSpline(const std::array<KnotVector, DIM> &knot_vector,
          std::array<int, DIM> degree,
          const std::vector<ControlPoint> &control_points)
      : dim(control_points[0].GetDimension()) {
    for (int i = 0; i < DIM, ++i) {
      parameter_space_[i] = parameter_space_(knot_vector[i], degree[i]);
    }
    for (auto &&cp : control_points) {
      for (int i = 0; i < dim; ++i) {
        control_points_.emplace_back(cp.GetValue(i));
      }
    }
  }

  std::vector<double> Evaluate(std::array<double, DIM> param_coord, const std::vector<int> &dimensions) const {
    auto basis_function_values = parameter_space_.EvaluateAllNonZeroBasisFunctions(param_coord);
    std::vector<double> evaluated_point(dimensions.size(), 0);
    for (int i = 0; i < dimensions.size(); ++i) {
      evaluated_point[i] =
          ComputeWeightedSum(basis_function_values, ExtractControlPointValues(param_coord, dimensions[i]));
    }
    return evaluated_point;
  }

  std::vector<double> EvaluateDerivative(double param_coord, const std::vector<int> &dimensions, int derivative) const {
    auto basis_function_values = parameter_space_.EvaluateAllNonZeroBasisFunctionDerivatives(param_coord, derivative);
    std::vector<double> evaluated_point(dimensions.size(), 0);
    for (int i = 0; i < dimensions.size(); ++i) {
      evaluated_point[i] =
          ComputeWeightedSum(basis_function_values, ExtractControlPointValues(param_coord, dimensions[i]));
    }
    return evaluated_point;
  }

  int GetDegree() const {
    return parameter_space_.degree();
  }

  KnotVector GetKnotVector() const {
    return parameter_space_.knot_vector();
  }

 private:
  std::vector<double> ExtractControlPointValues(double param_coord, int dimension) const {
    std::vector<double> control_point_values(static_cast<uint64_t>(GetDegree() + 1), 0.0);
    auto control_points =
        control_points_.begin() + (GetKnotVector().GetKnotSpan(param_coord) - GetDegree())*dim + dimension;
    for (int i = 0; i < GetDegree() + 1; ++i) {
      control_point_values[i] = *control_points;
      control_points += dim;
    }
    return control_point_values;
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

  std::vector<double> EvaluateAllNonZeroBasisFunctions(std::array<double, DIM> param_coord) const {
    std::array<std::vector<std::unique_ptr<BasisFunction>>::const_iterator, DIM> first_non_zero;
    for (int i = 0; i < DIM; ++i) {
      first_non_zero[i] = parameter_space_[i].GetFirstNonZeroKnot(param_coord[i]);
    }
    std::array<int, DIM> lastKnotOffset;
    for (int i = 0; i < DIM; ++i) {
      lastKnotOffset[i] = parameter_space_[i].degree_;
    }
    MultiIndexHandler MultiIndexHandler(lastKnotOffset);
    int M = 1;
    for (int i = 0; i < DIM; ++i) {
      M *= tensor[i].size();
    }
    std::vector<double> vector(M, 1);
    for (int i = 0; i < M; ++i) {
      for (int j = 0; j < DIM; ++j) {
        vector[i] *= (first_non_zero[j] + MultiIndexHandler[j])->Evaluate(param_coord[j]);
        MultiIndexHandler++;
      }
    }
    return vector;
  }

  std::array<ParameterSpace, DIM> parameter_space_;
  std::vector<double> control_points_;
  int dim;
};

#endif  // SRC_B_SPLINE_H_
