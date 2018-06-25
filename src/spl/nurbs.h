/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_SPL_NURBS_H
#define SRC_SPL_NURBS_H

#include <algorithm>
#include <array>
#include <functional>
#include <utility>
#include <vector>

#include "b_spline.h"
#include "spline.h"
#include "weighted_physical_space.h"

namespace spl {
template<int DIM>
class NURBS : public Spline<DIM> {
 public:
  NURBS(std::shared_ptr<std::array<baf::KnotVector, DIM >> knot_vector,
        std::array<int, DIM> degree,
        const std::vector<baf::ControlPoint> &control_points, std::vector<double> weights) : Spline<DIM>(knot_vector,
                                                                                                         degree) {
    std::array<int, DIM> number_of_points;
    for (int i = 0; i < DIM; ++i) {
      number_of_points[i] = (*knot_vector)[i].GetNumberOfKnots() - degree[i] - 1;
    }
    physical_space_ = WeightedPhysicalSpace<DIM>(control_points, weights, number_of_points);
  }

  NURBS(ParameterSpace<DIM> parameter_space, WeightedPhysicalSpace<DIM> physical_space) : Spline<DIM>(std::move(
      parameter_space)), physical_space_(physical_space) {}

 private:
  double GetEvaluatedControlPoint(std::array<ParamCoord, DIM> param_coord,
                                  std::array<int, DIM> indices,
                                  int dimension) const override {
    return this->parameter_space_.GetBasisFunctions(indices, param_coord)
        * physical_space_.GetHomogenousControlPoint(indices).GetValue(dimension)
        / GetEvaluatedDerivativeWeightSum(param_coord, std::array<int, DIM>{0});
  }

  double GetEvaluatedDerivativeControlPoint(std::array<ParamCoord, DIM> param_coord,
                                            std::array<int, DIM> derivative,
                                            std::array<int, DIM> indices,
                                            int dimension) const override {
    return GetRationalBasisFunctionDerivative(param_coord, derivative, indices, dimension)
        * physical_space_.GetControlPoint(indices).GetValue(dimension);
  }

  double GetRationalBasisFunctionDerivative(std::array<ParamCoord, DIM> param_coord,
                                            std::array<int, DIM> derivative,
                                            std::array<int, DIM> indices,
                                            int dimension) const {
    if (derivative == std::array<int, DIM>{0}) {
      return this->parameter_space_.GetBasisFunctions(indices, param_coord)
          * physical_space_.GetWeight(indices) / GetEvaluatedDerivativeWeightSum(param_coord, std::array<int, DIM>{0});
    }
    return (GetEvaluatedDerivativeWeight(param_coord, derivative, indices)
        - GetDerivativesSum(param_coord, derivative, indices, dimension))
        / GetEvaluatedDerivativeWeightSum(param_coord, std::array<int, DIM>{0});
  }

  double GetDerivativesSum(std::array<ParamCoord, DIM> param_coord,
                           std::array<int, DIM> derivative,
                           std::array<int, DIM> indices,
                           int dimension) const {
    util::MultiIndexHandler<DIM> derivativeHandler(GetDerivativeHandler(derivative));
    derivativeHandler++;
    double sum = 0;
    for (int i = 1; i < derivativeHandler.Get1DLength(); i++, derivativeHandler++) {
      sum += binomialCoefficient(derivative, derivativeHandler.GetIndices())
          * GetRationalBasisFunctionDerivative(param_coord,
                                               derivativeHandler.GetDifferenceIndices(),
                                               indices,
                                               dimension)
          * GetEvaluatedDerivativeWeightSum(param_coord, derivativeHandler.GetIndices());
    }
    return sum;
  }

  double GetEvaluatedDerivativeWeightSum(std::array<ParamCoord, DIM> param_coord,
                                         std::array<int, DIM> derivative) const {
    auto first_non_zero = this->GetArrayOfFirstNonZeroBasisFunctions(param_coord);
    util::MultiIndexHandler<DIM> basisFunctionHandler(this->GetNumberOfBasisFunctionsToEvaluate());
    double sum = 0;
    for (int i = 0; i < basisFunctionHandler.Get1DLength(); ++i, basisFunctionHandler++) {
      auto indices = basisFunctionHandler.GetIndices();
      std::transform(indices.begin(), indices.end(), first_non_zero.begin(), indices.begin(), std::plus<double>());
      sum += GetEvaluatedDerivativeWeight(param_coord, derivative, indices);
    }
    return sum;
  }

  double GetEvaluatedDerivativeWeight(std::array<ParamCoord, DIM> param_coord,
                                      std::array<int, DIM> derivative,
                                      std::array<int, DIM> indices) const {
    return this->parameter_space_.GetBasisFunctionDerivatives(indices, param_coord, derivative)
        * physical_space_.GetWeight(indices);
  }

  util::MultiIndexHandler<DIM> GetDerivativeHandler(const std::array<int, DIM> &derivative) const {
    std::array<int, DIM> derivative_length;
    for (int i = 0; i < DIM; ++i) {
      derivative_length[i] = derivative[i] + 1;
    }
    return util::MultiIndexHandler<DIM>(derivative_length);
  }

  int binomialCoefficient(int number, int subset) const {
    if (subset == 0 || subset == number)
      return 1;
    return binomialCoefficient(number - 1, subset - 1) + binomialCoefficient(number - 1, subset);
  }

  int binomialCoefficient(std::array<int, DIM> numbers, std::array<int, DIM> subsets) const {
    int bc = 1;
    for (int i = 0; i < DIM; ++i) {
      bc *= binomialCoefficient(numbers[i], subsets[i]);
    }
    return bc;
  }

  WeightedPhysicalSpace<DIM> physical_space_;
};
}  // namespace spl

#endif  // SRC_SPL__NURBS_H
