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

template <int DIM>
class BSpline {
 public:
  BSpline(const std::array<KnotVector, DIM> &knot_vector,
          std::array<int, DIM> degree,
          const std::vector<ControlPoint> &control_points)
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
                                         int derivative) const {
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

  KnotVector GetKnotVector(int i) const {
    return parameter_space_[i].knot_vector();
  }

  std::vector<Element> GetElementList() const;
  std::vector<std::vector<double>> EvaluateAllElementNonZeroBasisFunctions(int element_number,
                                                                           const IntegrationRule<1> &rule) const;
  std::vector<std::vector<double>> EvaluateAllElementNonZeroBasisFunctionDerivatives(
      int element_number,
      const IntegrationRule<1> &rule) const;

  double JacobianDeterminant(int element_number, int integration_point, const IntegrationRule<1> &rule) const;

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
      total_length[i + 1] = parameter_space_[i].knot_vector().Size() - parameter_space_[i].degree() - 1;
      current[i + 1] = start[i + 1];
      M *= (last[i + 1] - start[i + 1]);
    }
    MultiIndexHandler<DIM + 1> multiIndexHandler(total_length);
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

  std::vector<std::vector<double>> TransformToPhysicalSpace(std::vector<std::vector<double>> values,
                                                            int element_number,
                                                            const IntegrationRule<1> &rule) const;
  double TransformToParameterSpace(double upper, double lower, double point) const;

  std::vector<double> EvaluateAllNonZeroBasisFunctions(std::array<double, DIM> param_coord) const {
    std::array<std::vector<std::unique_ptr<BasisFunction>>::const_iterator, DIM> first_non_zero;
    for (int i = 0; i < DIM; ++i) {
      first_non_zero[i] = parameter_space_[i].GetFirstNonZeroKnot(param_coord[i]);
    }
    std::array<int, DIM> total_length;
    for (int i = 0; i < DIM; ++i) {
      total_length[i] = parameter_space_[i].degree() + 1;
    }
    MultiIndexHandler<DIM> multiIndexHandler(total_length);
    int M = 1;
    for (int i = 0; i < DIM; ++i) {
      M *= total_length[i];
    }
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
                                                                 int derivative) const {
    std::array<std::vector<std::unique_ptr<BasisFunction>>::const_iterator, DIM> first_non_zero;
    for (int i = 0; i < DIM; ++i) {
      first_non_zero[i] = parameter_space_[i].GetFirstNonZeroKnot(param_coord[i]);
    }
    std::array<int, DIM> total_length;
    for (int i = 0; i < DIM; ++i) {
      total_length[i] = parameter_space_[i].degree() + 1;
    }
    MultiIndexHandler<DIM> multiIndexHandler(total_length);
    int M = 1;
    for (int i = 0; i < DIM; ++i) {
      M *= total_length[i];
    }
    std::vector<double> vector(M, 1);
    for (int i = 0; i < M; ++i) {
      for (int j = 0; j < DIM; ++j) {
        vector[i] *= (*(first_non_zero[j] + multiIndexHandler[j]))->EvaluateDerivative(derivative, param_coord[j]);
      }
      multiIndexHandler++;
    }
    return vector;
  }

  std::array<ParameterSpace, DIM> parameter_space_;
  std::vector<double> control_points_;
  int dim;
};

#endif  // SRC_B_SPLINE_H_
