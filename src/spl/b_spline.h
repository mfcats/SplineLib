/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_SPL_B_SPLINE_H_
#define SRC_SPL_B_SPLINE_H_

#include <algorithm>
#include <array>
#include <functional>
#include <vector>

#include "spline.h"

namespace spl {
template<int DIM>
class BSpline : public Spline<DIM> {
 public:
  BSpline(std::shared_ptr<std::array<baf::KnotVector, DIM>> knot_vector,
          std::array<int, DIM> degree,
          const std::vector<baf::ControlPoint> &control_points) : Spline<DIM>(knot_vector, degree) {
    std::array<int, DIM> number_of_points;
    for (int i = 0; i < DIM; ++i) {
      number_of_points[i] = (*knot_vector)[i].GetNumberOfKnots() - degree[i] - 1;
    }
    physical_space_ = PhysicalSpace<DIM>(control_points, number_of_points);
  }

  BSpline(ParameterSpace<DIM> parameter_space, PhysicalSpace<DIM> physical_space)
      : Spline<DIM>(std::move(parameter_space)), physical_space_(physical_space) {
  }

  std::vector<double> EvaluateDerivative(std::array<ParamCoord, DIM> param_coord,
                                         const std::vector<int> &dimensions,
                                         std::array<int, DIM> derivative) const override {
    this->ThrowIfParametricCoordinateOutsideKnotVectorRange(param_coord);
    auto basis_function_values = EvaluateAllNonZeroBasisFunctionDerivatives(param_coord, derivative);
    std::vector<double> evaluated_point(dimensions.size(), 0);
    for (int i = 0; i < dimensions.size(); ++i) {
      evaluated_point[i] =
          this->ComputeWeightedSum(basis_function_values, this->ExtractControlPointValues(param_coord, dimensions[i]));
    }
    return evaluated_point;
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
      start[i + 1] = this->GetKnotVector(i).GetKnotSpan(param_coord[i]) - this->GetDegree(i);
      last[i + 1] = start[i + 1] + this->parameter_space_.degree(i) + 1;
      total_length[i + 1] =
          this->parameter_space_.knot_vector(i).GetNumberOfKnots() - this->parameter_space_.degree(i) - 1;
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

  std::vector<double> EvaluateAllNonZeroBasisFunctionDerivatives(std::array<ParamCoord, DIM> param_coord,
                                                                 std::array<int, DIM> derivative) const {
    auto first_non_zero = this->CreateArrayFirstNonZeroBasisFunction(param_coord);
    auto total_length = this->GetNumberOfBasisFunctionsToEvaluate();

    util::MultiIndexHandler<DIM> multiIndexHandler(total_length);
    auto M = multiIndexHandler.Get1DLength();

    std::vector<double> vector(M, 1);
    for (int i = 0; i < M; ++i) {
      for (int j = 0; j < DIM; ++j) {
        vector[i] *= (*(first_non_zero[j] + multiIndexHandler[j]))->EvaluateDerivative(param_coord[j], derivative[j]);
      }
      multiIndexHandler++;
    }
    return vector;
  }

 private:
  double ComputeWeightedSum(const std::vector<double> &basis_function_values,
                            std::vector<double> control_point_values) const {
    std::transform(basis_function_values.begin(),
                   basis_function_values.end(),
                   control_point_values.begin(),
                   control_point_values.begin(),
                   std::multiplies<double>());
    return std::accumulate(control_point_values.begin(), control_point_values.end(), 0.0, std::plus<double>());
  }

  std::array<std::vector<std::shared_ptr<baf::BasisFunction>>::const_iterator, DIM>
  CreateArrayFirstNonZeroBasisFunction(std::array<ParamCoord, DIM> param_coord) const {
    std::array<std::vector<std::shared_ptr<baf::BasisFunction>>::const_iterator, DIM> first_non_zero;
    for (int i = 0; i < DIM; ++i) {
      first_non_zero[i] = this->parameter_space_.GetFirstNonZeroKnot(i, param_coord[i]);
    }
    return first_non_zero;
  }

  double GetEvaluatedControlPoint(std::array<ParamCoord, DIM> param_coord,
                                  std::array<int, DIM> indices,
                                  int dimension) const {
    return this->parameter_space_.GetBasisFunctions(indices, param_coord)
        * physical_space_.GetControlPoint(indices).GetValue(dimension);
  }

  PhysicalSpace<DIM> physical_space_;
};
}  // namespace spl

#endif  // SRC_SPL_B_SPLINE_H_
