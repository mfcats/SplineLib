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
          const std::vector<baf::ControlPoint> &control_points) : Spline<DIM>(knot_vector, degree, control_points) {}

  BSpline(ParameterSpace<DIM> parameter_spaces, PhysicalSpace<DIM> physical_space)
      : Spline<DIM>(std::move(parameter_spaces), physical_space) {
  }

  std::vector<double> Evaluate(std::array<ParamCoord, DIM> param_coord,
                               const std::vector<int> &dimensions) const override {
    this->ThrowIfParametricCoordinateOutsideKnotVectorRange(param_coord);

    std::array<int, DIM> first_non_zero;
    for (int i = 0; i < DIM; ++i) {
      first_non_zero[i] =
          this->parameter_space_.knot_vector(i).GetKnotSpan(param_coord[i]) - this->parameter_space_.degree(i);
    }

    auto total_length = this->ArrayTotalLength();
    auto M = this->MultiIndexHandlerShort();

    util::MultiIndexHandler<DIM> multiIndexHandler(total_length);
    std::vector<double> vector(dimensions.size(), 0);
    for (int i = 0; i < M; ++i) {
      for (int j = 0; j < dimensions.size(); ++j) {
        auto a = multiIndexHandler.GetIndices();
        std::transform(a.begin(), a.end(), first_non_zero.begin(), a.begin(), std::plus<double>());
        vector[j] += this->parameter_space_.GetBasisFunctions(a, param_coord)
            * this->physical_space_->GetControlPoint(a).GetValue(dimensions[j]);
      }
      multiIndexHandler++;
    }
    return vector;
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

  std::vector<double> EvaluateAllNonZeroBasisFunctionDerivatives(std::array<ParamCoord, DIM> param_coord,
                                                                 std::array<int, DIM> derivative) const {
    auto first_non_zero = this->CreateArrayFirstNonZeroBasisFunction(param_coord);
    auto total_length = this->ArrayTotalLength();
    auto M = this->MultiIndexHandlerShort();

    util::MultiIndexHandler<DIM> multiIndexHandler(total_length);

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
  std::vector<double> EvaluateAllNonZeroBasisFunctions(std::array<ParamCoord, DIM> param_coord) const override {
    auto first_non_zero = this->CreateArrayFirstNonZeroBasisFunction(param_coord);
    auto total_length = this->ArrayTotalLength();
    auto M = this->MultiIndexHandlerShort();

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
};
}  // namespace spl

#endif  // SRC_SPL_B_SPLINE_H_
