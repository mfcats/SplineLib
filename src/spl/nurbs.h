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

#include <array>
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

  std::vector<double> EvaluateDerivative(std::array<ParamCoord, DIM> param_coord,
                                         const std::vector<int> &dimensions,
                                         std::array<int, DIM> derivative) const override {
    this->ThrowIfParametricCoordinateOutsideKnotVectorRange(param_coord);
    if (derivative == std::array<int, DIM>{0}) {
      return this->Evaluate(param_coord, dimensions);
    }
    std::vector<double> evaluated_point(dimensions.size(), 0);
    for (int i = 0; i < dimensions.size(); ++i) {
      evaluated_point[i] = (GetHomogenousDerivative(param_coord, dimensions[i], derivative)
          - GetSum(param_coord, derivative, dimensions[i]))
          / GetWeightDerivative(param_coord, std::array<int, DIM>{0});
    }
    return evaluated_point;
  }

 private:
  double GetEvaluatedControlPoint(std::array<ParamCoord, DIM> param_coord,
                                  std::array<int, DIM> indices,
                                  int dimension) const override {
    return this->parameter_space_.GetBasisFunctions(indices, param_coord)
        * physical_space_.GetControlPoint(indices).GetValue(dimension) / GetSum(param_coord);
  }

  util::MultiIndexHandler<DIM> GetDerivativeHandler(const std::array<int, DIM> &derivative) const {
    std::array<int, DIM> derivative_length;
    for (int i = 0; i < DIM; ++i) {
      derivative_length[i] = derivative[i] + 1;
    }
    return util::MultiIndexHandler<DIM>(derivative_length);
  }

  double GetSum(const std::array<ParamCoord, DIM> &param_coord,
                const std::array<int, DIM> &derivative,
                int dimension) const {
    double sum = 0;
    util::MultiIndexHandler<DIM> derivativeHandler = GetDerivativeHandler(derivative);
    derivativeHandler++;
    for (int i = 0; i < derivativeHandler.Get1DLength() - 1; i++, derivativeHandler++) {
      sum += binomialCoefficient(derivative, derivativeHandler.GetIndices())
          * GetWeightDerivative(param_coord, derivativeHandler.GetIndices())
          * EvaluateDerivative(param_coord, {dimension}, derivativeHandler.GetDifferenceIndices())[0];
    }
    return sum;
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
    for (int control_point = 0; control_point < physical_space_.GetNumberOfControlPoints(); ++control_point) {
      auto a = physical_space_.GetWeight(control_point);
      weights.emplace_back(baf::ControlPoint({physical_space_.GetWeight(control_point)}));
    }
    return weights;
  }

  std::vector<baf::ControlPoint> GetHomogenousControlPoints() const {
    std::vector<baf::ControlPoint> homogenousPoints;
    for (int point = 0; point < physical_space_.GetNumberOfControlPoints(); ++point) {
      std::vector<double> homogenousCoordinates;
      for (int coordinate = 0; coordinate < physical_space_.GetDimension(); ++coordinate) {
        homogenousCoordinates.emplace_back(
            physical_space_.GetPoint(point * physical_space_.GetDimension() + coordinate)
                * physical_space_.GetWeight(point));
      }
      homogenousPoints.emplace_back(baf::ControlPoint(homogenousCoordinates));
    }
    return homogenousPoints;
  }

  std::unique_ptr<spl::BSpline<DIM>> GetBSpline(std::vector<baf::ControlPoint> controlPoints) const {
    std::shared_ptr<std::array<baf::KnotVector, DIM>>
        knot_vector_ptr = std::make_shared<std::array<baf::KnotVector, DIM>>(GetKnotVectors());
    return std::make_unique<spl::BSpline<DIM>>(knot_vector_ptr, GetDegrees(), controlPoints);
  }

  double GetSum(std::array<ParamCoord, DIM> param_coord) const {
    return GetBSpline(GetWeightsAsControlPoints())->Evaluate(param_coord, {0})[0];
  }

  double GetWeightDerivative(std::array<ParamCoord, DIM> param_coord, std::array<int, DIM> derivative) const {
    return GetBSpline(GetWeightsAsControlPoints())->EvaluateDerivative(param_coord, {0}, derivative)[0];
  }

  double GetHomogenousDerivative(std::array<ParamCoord, DIM> param_coord,
                                 int dimension,
                                 std::array<int, DIM> derivative) const {
    return GetBSpline(GetHomogenousControlPoints())->EvaluateDerivative(param_coord, {dimension}, derivative)[0];
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
