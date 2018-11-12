/* Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.
*/

#ifndef SRC_SPL_NURBS_H_
#define SRC_SPL_NURBS_H_

#include <algorithm>
#include <array>
#include <functional>
#include <utility>
#include <vector>
#include <iostream>

#include "b_spline.h"
#include "spline.h"
#include "nurbs_generator.h"
#include "weighted_physical_space.h"

namespace spl {
template<int DIM>
class NURBS : public Spline<DIM> {
 public:
  NURBS(std::array<std::shared_ptr<baf::KnotVector>, DIM> knot_vector,
        std::array<Degree, DIM> degree,
        const std::vector<baf::ControlPoint> &control_points, std::vector<double> weights) : Spline<DIM>(knot_vector,
                                                                                                         degree) {
    std::array<int, DIM> number_of_points;
    for (int i = 0; i < DIM; ++i) {
      number_of_points[i] = knot_vector[i]->GetNumberOfKnots() - degree[i].get() - 1;
    }
    physical_space_ = std::make_shared<WeightedPhysicalSpace<DIM>>(WeightedPhysicalSpace<DIM>(control_points,
                                                                                              weights,
                                                                                              number_of_points));
  }

  virtual ~NURBS() = default;

  NURBS(std::shared_ptr<ParameterSpace<DIM>> parameter_space,
        std::shared_ptr<WeightedPhysicalSpace<DIM>> physical_space) :
      physical_space_(physical_space) {
    this->parameter_space_ = parameter_space;
  }

  int GetNumberOfControlPoints() override {
    return physical_space_->GetNumberOfControlPoints();
  }

  std::array<int, DIM> GetPointsPerDirection() override {
    return physical_space_->GetNumberOfPointsInEachDirection();
  }

  int GetDimension() override {
    return physical_space_->GetDimension();
  }

  double GetControlPoint(std::array<int, DIM> indices, int dimension) override {
    return physical_space_->GetControlPoint(indices).GetValue(dimension);
  }

  double GetWeight(std::array<int, DIM> indices) {
    return physical_space_->GetWeight(indices);
  }

  explicit NURBS(NURBSGenerator<DIM> nurbs_generator) : Spline<DIM>(nurbs_generator.GetParameterSpace()) {
    physical_space_ = nurbs_generator.GetWeightedPhysicalSpace();
  }

  std::shared_ptr<spl::PhysicalSpace<DIM>> GetPhysicalSpace() const override {
    return physical_space_;
  }

  void AdjustControlPoints(std::vector<double> scaling, int first, int last, int dimension) override {
    std::array<int, DIM> number_of_points = physical_space_->GetNumberOfPointsInEachDirection();
    for (auto &number : number_of_points) {
      --number;
    }
    util::MultiIndexHandler<DIM> point_handler(physical_space_->GetNumberOfPointsInEachDirection());
    point_handler.SetIndices(number_of_points);
    for (int i = point_handler.Get1DLength() - 1; i >= 0; --i, --point_handler) {
      auto current_point = point_handler.GetIndices()[dimension];
      if (current_point <= last && point_handler.GetIndices()[dimension] >= first) {
        std::array<int, DIM> indices0 = point_handler.GetIndices(), indices1 = indices0;
        --indices1[dimension];
        baf::ControlPoint cp0 = physical_space_->GetHomogenousControlPoint(indices0);
        baf::ControlPoint cp1 = physical_space_->GetHomogenousControlPoint(indices1);
        double weight0 = physical_space_->GetWeight(indices0);
        double weight1 = physical_space_->GetWeight(indices1);
        double new_weight = scaling[i - first] * weight0 + (1 - scaling[i - first]) * weight1;
        std::vector<double> coordinates;
        for (int j = 0; j < cp0.GetDimension(); ++j) {
          coordinates.push_back((scaling[current_point - first] * cp0.GetValue(j)
              + (1 - scaling[current_point - first]) * cp1.GetValue(j)) / new_weight);
        }
        baf::ControlPoint new_cp(coordinates);
        current_point != last ? (physical_space_->SetControlPoint(indices0, new_cp, dimension),
            physical_space_->SetWeight(indices0, new_weight))
                              : (physical_space_->InsertControlPoint(indices0, new_cp),
            physical_space_->InsertWeight(indices0, new_weight));
      }
    }
  }

 private:
  double GetEvaluatedControlPoint(std::array<ParamCoord, DIM> param_coord,
                                  std::array<int, DIM> indices,
                                  int dimension) const override {
    return this->parameter_space_->GetBasisFunctions(indices, param_coord)
        * physical_space_->GetHomogenousControlPoint(indices).GetValue(dimension)
        / GetEvaluatedWeightSum(param_coord);
  }

  double GetEvaluatedDerivativeControlPoint(std::array<ParamCoord, DIM> param_coord,
                                            std::array<int, DIM> derivative,
                                            std::array<int, DIM> indices,
                                            int dimension) const override {
    return GetRationalBasisFunctionDerivative(param_coord, derivative, indices, dimension)
        * physical_space_->GetControlPoint(indices).GetValue(dimension);
  }

  double GetRationalBasisFunctionDerivative(std::array<ParamCoord, DIM> param_coord,
                                            std::array<int, DIM> derivative,
                                            std::array<int, DIM> indices,
                                            int dimension) const {
    if (derivative == std::array<int, DIM>{0}) {
      return this->parameter_space_->GetBasisFunctions(indices, param_coord)
          * physical_space_->GetWeight(indices) / GetEvaluatedWeightSum(param_coord);
    }
    return (GetEvaluatedDerivativeWeight(param_coord, derivative, indices)
        - GetDerivativesSum(param_coord, derivative, indices, dimension))
        / GetEvaluatedWeightSum(param_coord);
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

  double GetEvaluatedWeightSum(std::array<ParamCoord, DIM> param_coord) const {
    auto first_non_zero = this->GetArrayOfFirstNonZeroBasisFunctions(param_coord);
    util::MultiIndexHandler<DIM> basisFunctionHandler(this->GetNumberOfBasisFunctionsToEvaluate());
    double sum = 0;
    for (int i = 0; i < basisFunctionHandler.Get1DLength(); ++i, basisFunctionHandler++) {
      auto indices = basisFunctionHandler.GetIndices();
      std::transform(indices.begin(), indices.end(), first_non_zero.begin(), indices.begin(), std::plus<double>());
      sum += GetEvaluatedWeight(param_coord, indices);
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

  double GetEvaluatedWeight(std::array<ParamCoord, DIM> param_coord, std::array<int, DIM> indices) const {
    return this->parameter_space_->GetBasisFunctions(indices, param_coord) * physical_space_->GetWeight(indices);
  }

  double GetEvaluatedDerivativeWeight(std::array<ParamCoord, DIM> param_coord,
                                      std::array<int, DIM> derivative,
                                      std::array<int, DIM> indices) const {
    return this->parameter_space_->GetBasisFunctionDerivatives(indices, param_coord, derivative)
        * physical_space_->GetWeight(indices);
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

  std::shared_ptr<WeightedPhysicalSpace<DIM>> physical_space_;
};
}  // namespace spl

#endif  // SRC_SPL_NURBS_H_
