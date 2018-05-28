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
    if (DIM > 2) {
      throw std::runtime_error("Evaluation of NURBS derivatives only implemented for NURBS of dimensions 1 and 2");
    }
    if (derivative == std::array<int, DIM>{0}) {
      return this->Evaluate(param_coord, dimensions);
    }
    std::array<int, DIM> derivative_length;
    for (int i = 0; i < DIM; ++i) {
      derivative_length[i] = derivative[i] + 1;
    }
    util::MultiIndexHandler<DIM> derivativeHandler(derivative_length);
    std::vector<double> evaluated_point(dimensions.size(), 0);
    for (int i = 0; i < dimensions.size(); ++i) {
      double sum = 0;
      for (int j = 1; j <= derivative[0]; ++j) {
        derivativeHandler++;
        sum += binomialCoefficient(derivative[0], j) * GetWeightDerivative(param_coord, derivativeHandler.GetIndices())
            * EvaluateDerivative(param_coord, {dimensions[i]}, derivativeHandler.GetDifferenceIndices())[0];
      }
      if (DIM > 1) {
        std::array<int, DIM> zeros = {0};
        derivativeHandler.SetIndices(zeros);
        for (int j = 1; j <= derivative[1]; ++j) {
          derivativeHandler + derivative[0] + 1;
          sum +=
              binomialCoefficient(derivative[1], j) * GetWeightDerivative(param_coord, derivativeHandler.GetIndices())
                  * EvaluateDerivative(param_coord, {dimensions[i]}, derivativeHandler.GetDifferenceIndices())[0];
        }
        derivativeHandler.SetIndices(zeros);
        for (int k = 1; k <= derivative[0]; ++k) {
          derivativeHandler++;
          int sum2 = binomialCoefficient(derivative[0], k);
          for (int l = 1; l <= derivative[1]; ++l) {
            derivativeHandler + derivative[0] + 1;
            sum2 *=
                binomialCoefficient(derivative[1], l) * GetWeightDerivative(param_coord, derivativeHandler.GetIndices())
                    * EvaluateDerivative(param_coord, {dimensions[i]}, derivativeHandler.GetDifferenceIndices())[0];
          }
          if (derivative[1] > 0) {
            sum += sum2;
          }
        }
      }
      evaluated_point[i] = (GetHomogenousDerivative(param_coord, dimensions[i], derivative) - sum)
          / GetWeightDerivative(param_coord, std::array<int, DIM>{0});
    }
    return evaluated_point;
  }

 private:
  std::vector<double> EvaluateAllNonZeroBasisFunctions(std::array<double, DIM> param_coord) const override {
    auto first_non_zero = this->CreateArrayFirstNonZeroBasisFunction(param_coord);
    util::MultiIndexHandler<DIM> multiIndexHandler(this->ArrayTotalLength());
    std::vector<double> NonZeroBasisFunctions(this->MultiIndexHandlerShort(), 1);
    auto extractedWeights = GetBSpline(GetWeightsAsControlPoints())->ExtractControlPointValues(param_coord, 0);
    for (double &basis_function : NonZeroBasisFunctions) {
      for (int j = 0; j < DIM; ++j) {
        basis_function *= (*(first_non_zero[j] + multiIndexHandler[j]))->Evaluate(param_coord[j]);
      }
      basis_function *= extractedWeights[multiIndexHandler.Get1DIndex()] / GetSum(param_coord);
      multiIndexHandler++;
    }
    return NonZeroBasisFunctions;
  }

  std::array<baf::KnotVector, DIM> GetKnotVectors() const {
    std::array<baf::KnotVector, DIM> knot_vectors;
    for (int vector = 0; vector < DIM; ++vector) {
      knot_vectors[vector] = this->GetKnotVector(vector);
    }
    return knot_vectors;
  }

  std::array<int, DIM> GetDegrees() const {
    std::array<int, DIM> degrees;
    for (int degree = 0; degree < DIM; ++degree) {
      degrees[degree] = this->GetDegree(degree);
    }
    return degrees;
  }

  std::vector<baf::ControlPoint> GetWeightsAsControlPoints() const {
    std::vector<baf::ControlPoint> weights;
    for (int control_point = 0; control_point < weights_.size(); ++control_point) {
      weights.emplace_back(baf::ControlPoint({weights_[control_point]}));
    }
    return weights;
  }

  std::vector<baf::ControlPoint> GetHomogenousControlPoints() const {
    std::vector<baf::ControlPoint> homogenousPoints;
    for (int point = 0; point < this->control_points_.size(); ++point) {
      std::vector<double> homogenousCoordinates;
      for (int coordinate = 0; coordinate < this->dim; ++coordinate) {
        homogenousCoordinates.emplace_back(this->control_points_[point * this->dim + coordinate] * weights_[point]);
      }
      homogenousPoints.emplace_back(baf::ControlPoint(homogenousCoordinates));
    }
    return homogenousPoints;
  }

  std::unique_ptr<spl::BSpline<DIM>> GetBSpline(std::vector<baf::ControlPoint> controlPoints) const {
    return std::make_unique<spl::BSpline<DIM>>(GetKnotVectors(), GetDegrees(), controlPoints);
  }

  double GetSum(std::array<double, DIM> param_coord) const {
    return GetBSpline(GetWeightsAsControlPoints())->Evaluate(param_coord, {0})[0];
  }

  double GetWeightDerivative(std::array<double, DIM> param_coord, std::array<int, DIM> derivative) const {
    return GetBSpline(GetWeightsAsControlPoints())->EvaluateDerivative(param_coord, {0}, derivative)[0];
  }

  double GetHomogenousDerivative(std::array<double, DIM> param_coord,
                                 int dimension,
                                 std::array<int, DIM> derivative) const {
    return GetBSpline(GetHomogenousControlPoints())->EvaluateDerivative(param_coord, {dimension}, derivative)[0];
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
