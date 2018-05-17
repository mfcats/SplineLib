/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SPLINELIB_NURBS_H
#define SPLINELIB_NURBS_H

#include <array>
#include <utility>
#include <vector>

#include "b_spline.h"
#include "spline.h"

namespace spl {
template<int DIM>
class NURBS : public Spline<DIM> {
 public:
  NURBS(const std::array<baf::KnotVector, DIM> &knot_vector,
        std::array<int, DIM> degree,
        const std::vector<baf::ControlPoint> &control_points, std::vector<double> weights) : Spline<DIM>(
      knot_vector,
      degree,
      control_points), weights_(std::move(weights)) {}

  std::vector<double> EvaluateDerivative(std::array<double, DIM> param_coord,
                                         const std::vector<int> &dimensions,
                                         std::array<int, DIM> derivative) const override {
    if (derivative[0] == 0) {
      auto test = this->Evaluate(param_coord, dimensions);
      return this->Evaluate(param_coord, dimensions);
    }
    std::vector<double> evaluated_point(dimensions.size(), 0);
    for (int i = 0; i < dimensions.size(); ++i) {
      double sum = 0;
      for (int j = 1; j <= derivative[0]; ++j) {
        auto a = GetWeightDerivative(param_coord, j);
        auto b = EvaluateDerivative(param_coord, {dimensions[i]}, {derivative[0] - j})[0];
        sum += binomialCoefficient(derivative[0], j) * GetWeightDerivative(param_coord, j)
            * EvaluateDerivative(param_coord, {dimensions[i]}, {derivative[0] - j})[0];
      }
      auto a = GetHomogenousDerivative(param_coord[0], dimensions[i], derivative[0]);
      auto b = GetWeightDerivative(param_coord, 0);
      evaluated_point[i] =
          (GetHomogenousDerivative(param_coord[0], dimensions[i], derivative[0]) - sum)
              / GetWeightDerivative(param_coord, 0);
    }
    return evaluated_point;
  }

 private:
  std::vector<double> EvaluateAllNonZeroBasisFunctions(std::array<double, DIM> param_coord) const override {
    auto first_non_zero = this->CreateArrayFirstNonZeroBasisFunction(param_coord);
    util::MultiIndexHandler<DIM> multiIndexHandler(this->ArrayTotalLength());
    std::vector<double> NonZeroBasisFunctions(this->MultiIndexHandlerShort(), 1);
    for (double &basis_function : NonZeroBasisFunctions) {
      for (int j = 0; j < DIM; ++j) {
        basis_function *= (*(first_non_zero[j] + multiIndexHandler[j]))->Evaluate(param_coord[j])
            * weights_[this->GetKnotVector(0).GetKnotSpan(param_coord[0]) - this->GetDegree(0) + multiIndexHandler[j]];
      }
      basis_function /= GetSum(param_coord);
      multiIndexHandler++;
    }
    return NonZeroBasisFunctions;
  }

  double GetSum(std::array<double, DIM> param_coord) const {
    std::vector<baf::ControlPoint> weights;
    for (int control_point = 0; control_point < weights_.size(); ++control_point) {
      weights.emplace_back(baf::ControlPoint({weights_[control_point]}));
    }
    return spl::BSpline<DIM>(std::array<baf::KnotVector, DIM>{this->GetKnotVector(0)},
                             std::array<int, DIM>{this->GetDegree(0)},
                             weights).Evaluate(param_coord, {0})[0];
  }

  double GetWeightDerivative(std::array<double, DIM> param_coord, int derivative) const {
    std::vector<baf::ControlPoint> weights;
    for (int control_point = 0; control_point < weights_.size(); ++control_point) {
      weights.emplace_back(baf::ControlPoint({weights_[control_point]}));
    }
    return BSpline<1>(std::array<baf::KnotVector, 1>{this->GetKnotVector(0)},
                      std::array<int, 1>{this->GetDegree(0)},
                      weights).EvaluateDerivative({param_coord[0]}, {0}, {derivative})[0];
  }

  double GetHomogenousDerivative(double param_coord, int dimension, int derivative) const {
    std::vector<baf::ControlPoint> homogenousPoints;
    for (int point = 0; point < this->control_points_.size(); ++point) {
      std::vector<double> homogenousCoordinates;
      for (int coordinate = 0; coordinate < this->dim; ++coordinate) {
        homogenousCoordinates.emplace_back(this->control_points_[point * this->dim + coordinate] * weights_[point]);
      }
      homogenousPoints.emplace_back(baf::ControlPoint(homogenousCoordinates));
    }
    return BSpline<1>(std::array<baf::KnotVector, 1>{this->GetKnotVector(0)},
                      std::array<int, 1>{this->GetDegree(0)},
                      homogenousPoints).EvaluateDerivative({param_coord}, {dimension}, {derivative})[0];
  }

  int binomialCoefficient(int number, int subset) const {
    if (subset == 0 || subset == number)
      return 1;
    return binomialCoefficient(number - 1, subset - 1) + binomialCoefficient(number - 1, subset);
  }

  std::vector<double> weights_;
};
} //namespace spl

#endif //SPLINELIB_NURBS_H
