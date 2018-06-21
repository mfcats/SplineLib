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

    auto first_non_zero = this->GetArrayOfFirstNonZeroBasisFunctions(param_coord);
    util::MultiIndexHandler<DIM> basisFunctionHandler(this->GetNumberOfBasisFunctionsToEvaluate());
    auto numberOfSummands = basisFunctionHandler.Get1DLength();
    std::vector<double> vector(dimensions.size(), 0);

    for (int i = 0; i < numberOfSummands; ++i) {
      auto indices = basisFunctionHandler.GetIndices();
      std::transform(indices.begin(), indices.end(), first_non_zero.begin(), indices.begin(), std::plus<double>());
      for (int j = 0; j < dimensions.size(); ++j) {
        vector[j] += GetEvaluatedDerivativeControlPoint(param_coord, derivative, indices, dimensions[j]);
      }
      basisFunctionHandler++;
    }
    return vector;
  }

 private:
  double GetEvaluatedControlPoint(std::array<ParamCoord, DIM> param_coord,
                                  std::array<int, DIM> indices,
                                  int dimension) const override {
    return this->parameter_space_.GetBasisFunctions(indices, param_coord)
        * physical_space_.GetControlPoint(indices).GetValue(dimension);
  }

  double GetEvaluatedDerivativeControlPoint(std::array<ParamCoord, DIM> param_coord,
                                            std::array<int, DIM> derivative,
                                            std::array<int, DIM> indices,
                                            int dimension) const {
    return this->parameter_space_.GetBasisFunctionDerivatives(indices, param_coord, derivative)
        * physical_space_.GetControlPoint(indices).GetValue(dimension);
  }

  PhysicalSpace<DIM> physical_space_;
};
}  // namespace spl

#endif  // SRC_SPL_B_SPLINE_H_
