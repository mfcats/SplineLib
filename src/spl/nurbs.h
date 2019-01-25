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

#include "b_spline.h"
#include "spline.h"
#include "nurbs_generator.h"
#include "weighted_physical_space.h"

namespace spl {
template<int DIM>
class NURBS : public Spline<DIM> {
 public:
  NURBS(KnotVectors<DIM> knot_vector, std::array<Degree, DIM> degree,
        const std::vector<baf::ControlPoint> &control_points,
        std::vector<double> weights) : Spline<DIM>(knot_vector, degree) {
    std::array<int, DIM> number_of_points;
    for (int i = 0; i < DIM; ++i) {
      number_of_points[i] = knot_vector[i]->GetNumberOfKnots() - degree[i].get() - 1;
    }
    physical_space_ = std::make_shared<WeightedPhysicalSpace<DIM>>(
        WeightedPhysicalSpace<DIM>(control_points, weights, number_of_points));
  }

  explicit NURBS(NURBSGenerator<DIM> nurbs_generator) : Spline<DIM>(nurbs_generator.GetParameterSpace()) {
    physical_space_ = nurbs_generator.GetWeightedPhysicalSpace();
  }

  NURBS(const NURBS<DIM> &nurbs) : Spline<DIM>(nurbs) {
    WeightedPhysicalSpace<DIM> weighted_physical_space(*nurbs.physical_space_);
    physical_space_ = std::make_shared<WeightedPhysicalSpace<DIM>>(weighted_physical_space);
  }

  virtual ~NURBS() = default;

  int GetNumberOfControlPoints() override {
    return physical_space_->GetNumberOfControlPoints();
  }

  std::array<int, DIM> GetPointsPerDirection() override {
    return physical_space_->GetNumberOfPointsInEachDirection();
  }

  int GetDimension() const override {
    return physical_space_->GetDimension();
  }

  double GetControlPoint(std::array<int, DIM> indices, int dimension) override {
    return physical_space_->GetControlPoint(indices).GetValue(dimension);
  }

  baf::ControlPoint GetControlPoint(std::array<int, DIM> indices) {
    return physical_space_->GetControlPoint(indices);
  }

  double GetWeight(std::array<int, DIM> indices) {
    return physical_space_->GetWeight(indices);
  }

  std::shared_ptr<spl::PhysicalSpace<DIM>> GetPhysicalSpace() const override {
    return physical_space_;
  }

//  std::shared_ptr<spl::WeightedPhysicalSpace<DIM>> GetWeightedPhysicalSpace() const {
//    return physical_space_;
//  }

  void AdjustControlPoints(std::vector<double> scaling, int first, int last, int dimension) override {
    std::array<int, DIM> point_handler_length = physical_space_->GetNumberOfPointsInEachDirection();
    ++point_handler_length[dimension];
    util::MultiIndexHandler<DIM> point_handler(point_handler_length);
    std::array<int, DIM> maximum_point_index = physical_space_->GetMaximumPointIndexInEachDirection();
    ++maximum_point_index[dimension];
    point_handler.SetIndices(maximum_point_index);
    int new_points = physical_space_->GetNumberOfControlPoints() / maximum_point_index[dimension];
    physical_space_->AddControlPoints(new_points);
    physical_space_->AddWeights(new_points);
    for (int i = point_handler.Get1DLength() - 1; i >= 0; --i, --point_handler) {
      auto current_point = point_handler[dimension];
      std::array<int, DIM> indices = point_handler.GetIndices();
      baf::ControlPoint new_control_point = GetNewControlPoint(indices, dimension, scaling, current_point, first, last);
      double new_weight = GetNewWeight(indices, dimension, scaling, current_point, first, last);
      physical_space_->SetControlPoint(indices, new_control_point, dimension);
      physical_space_->SetWeight(indices, new_weight, dimension);
    }
    physical_space_->IncrementNumberOfPoints(dimension);
  }

  bool RemoveControlPoints(std::vector<double> scaling, int first, int last, int dimension, double tolerance) override {
    int off = first - 1, i = first, j = last;
    std::vector<double> temp_w = GetTemporaryNewWeights(scaling, first, last, off, i, j);
    std::vector<double> temp = GetTemporaryNewControlPoints(scaling, temp_w, first, last, off, i, j);
    i += (j - i) / 2, j -= (j - i) / 2;
    if (!IsKnotRemovable(scaling[i - off - 1], temp, temp_w, tolerance, i, j, off, dimension)) {
      return false;
    }
    SetNewControlPoints(temp, last, i - off, off, dimension);
    SetNewWeights(temp_w, last, i - off, off, dimension);
    std::array<int, DIM> maximum_point_index = physical_space_->GetMaximumPointIndexInEachDirection();
    --maximum_point_index[dimension];
    physical_space_->RemoveControlPoints(
        physical_space_->GetNumberOfControlPoints() / physical_space_->GetNumberOfPointsInEachDirection()[dimension]);
    physical_space_->RemoveWeights(
        physical_space_->GetNumberOfControlPoints() / physical_space_->GetNumberOfPointsInEachDirection()[dimension]);
    physical_space_->DecrementNumberOfPoints(dimension);
    return true;
  }

  std::array<std::shared_ptr<spl::NURBS<DIM>>, 2> SudivideSpline(ParamCoord param_coord, int dimension) {
    this->InsertKnot(param_coord, dimension,
                     this->GetDegree(dimension).get() + 1
                         - this->GetKnotVector(dimension)->GetMultiplicity(param_coord));
    std::array<KnotVectors<DIM>, 2>
        new_knot_vectors = this->parameter_space_->GetDividedKnotVectors(param_coord, dimension);
    std::array<Degree, DIM> degrees = this->parameter_space_->GetDegrees();
    std::array<std::shared_ptr<spl::NURBS<DIM>>, 2> subdivided_splines;
    int first = 0;
    for (int i = 0; i < 2; ++i) {
      int length = new_knot_vectors[i][dimension]->GetNumberOfKnots() - degrees[dimension].get() - 1;
      std::vector<baf::ControlPoint> points = physical_space_->GetDividedControlPoints(first, length, dimension);
      std::vector<double> weights = physical_space_->GetDividedWeights(first, length, dimension);
      spl::NURBS<DIM> spline(new_knot_vectors[i], degrees, points, weights);
      subdivided_splines[i] = std::make_shared<spl::NURBS<DIM>>(spline);
      first = length;
    }
    return subdivided_splines;
  }

 private:
  double GetEvaluatedControlPoint(std::array<ParamCoord, DIM> param_coord,
                                  std::array<int, DIM> indices,
                                  int dimension) const override {
    if (GetDimension() == dimension) {
      return this->parameter_space_->GetBasisFunctions(indices, param_coord)
          * physical_space_->GetHomogenousControlPoint(indices).GetValue(dimension);
    } else {
      return this->parameter_space_->GetBasisFunctions(indices, param_coord)
          * physical_space_->GetHomogenousControlPoint(indices).GetValue(dimension)
          / GetEvaluatedWeightSum(param_coord);
    }
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

  baf::ControlPoint GetNewControlPoint(std::array<int, DIM> indices, int dimension, std::vector<double> scaling,
                                       int current_point, int first, int last) {
    if (current_point > last) {
      --indices[dimension];
      return physical_space_->GetControlPoint(indices);
    } else if (current_point >= first) {
      std::array<int, DIM> lower_indices = indices;
      --lower_indices[dimension];
      baf::ControlPoint upper_control_point = physical_space_->GetHomogenousControlPoint(indices);
      baf::ControlPoint lower_control_point = physical_space_->GetHomogenousControlPoint(lower_indices);
      std::vector<double> coordinates;
      double new_weight = GetNewWeight(indices, dimension, scaling, current_point, first, last);
      for (int j = 0; j < upper_control_point.GetDimension(); ++j) {
        coordinates.push_back((scaling[current_point - first] * upper_control_point.GetValue(j)
            + (1 - scaling[current_point - first]) * lower_control_point.GetValue(j)) / new_weight);
      }
      return baf::ControlPoint(coordinates);
    } else {
      return physical_space_->GetControlPoint(indices);
    }
  }

  double GetNewWeight(std::array<int, DIM> indices, int dimension, std::vector<double> scaling,
                      int current_point, int first, int last) {
    if (current_point > last) {
      --indices[dimension];
      return physical_space_->GetWeight(indices);
    } else if (current_point >= first) {
      std::array<int, DIM> lower_indices = indices;
      --lower_indices[dimension];
      double upper_weight = physical_space_->GetWeight(indices);
      double lower_weight = physical_space_->GetWeight(lower_indices);
      return scaling[current_point - first] * upper_weight + (1 - scaling[current_point - first]) * lower_weight;

    } else {
      return physical_space_->GetWeight(indices);
    }
  }

  void SetNewControlPoints(std::vector<double> temp, int last, int ii, int off, int dimension) {
    std::array<int, DIM> point_handler_length = physical_space_->GetNumberOfPointsInEachDirection();
    util::MultiIndexHandler<DIM> point_handler(point_handler_length);
    std::vector<double> coordinates(GetDimension(), 0);
    for (int m = 0; m < point_handler.Get1DLength(); ++m, ++point_handler) {
      int k = point_handler[dimension];
      if (k - off >= 1 && k - off != ii && k < last + 2) {
        int index = (point_handler.ExtractDimension(dimension) * (last - off + 2) + k - off) * GetDimension();
        for (int l = 0; l < GetDimension(); ++l) {
          coordinates[l] = temp[index + l];
        }
        auto indices = point_handler.GetIndices();
        indices[dimension] = k - off < ii ? k : k - 1;
        baf::ControlPoint cp(coordinates);
        physical_space_->SetControlPoint2(indices, cp, dimension);
      }
      if ((k <= off && k - off < 1)
          || (k >= last + 1 && k < physical_space_->GetNumberOfPointsInEachDirection()[dimension])) {
        auto indices = point_handler.GetIndices();
        indices[dimension] = k <= off ? k : k - 1;
        physical_space_->SetControlPoint2(indices,
                                          physical_space_->GetControlPoint(point_handler.GetIndices()),
                                          dimension);
      }
    }
  }

  void SetNewWeights(std::vector<double> temp, int last, int ii, int off, int dimension) {
    int index = 0;
    for (int k = 1; k < static_cast<int>(temp.size()); ++k) {
      if (k != ii) {
        index = k < ii ? k + off : k + off - 1;
        physical_space_->SetWeight2({index}, temp[k], dimension);
      }
    }
    for (int k = last + 1; k < physical_space_->GetNumberOfControlPoints(); ++k) {
      index = k - 1;
      physical_space_->SetWeight2({index}, physical_space_->GetWeights()[k], dimension);
    }
    auto weights = physical_space_->GetWeights();
  }

  std::vector<double> GetTemporaryNewControlPoints(std::vector<double> scaling, std::vector<double> temp_w, int first,
                                                   int last, int off, int i, int j) const {
    std::vector<double> temp((last + 2 - off) * GetDimension(), 0);
    for (int k = 0; k < GetDimension(); ++k) {
      temp[0 + k] = physical_space_->GetControlPoint({first - 1}).GetValue(k);
      temp[(last + 1 - off) * GetDimension() + k] =
          physical_space_->GetControlPoint({last + 1}).GetValue(k);
    }
    while (j - i > 0) {
      double alfi = scaling[i - first];
      double alfj = scaling[j - first];
      for (int k = 0; k < GetDimension(); ++k) {
        temp[(i - off) * GetDimension() + k] = (physical_space_->GetHomogenousControlPoint({i}).GetValue(k)
            - (1 - alfi) * temp[(i - off - 1) * GetDimension() + k]) / alfi / temp_w[i - off];
        temp[(j - off) * GetDimension() + k] =
            (physical_space_->GetHomogenousControlPoint({j}).GetValue(k)
                - alfj * temp[(j - off + 1) * GetDimension() + k])
                / (1 - alfj) / temp_w[j - off];
      }
      ++i, --j;
    }
    return temp;
  }

  std::vector<double>
  GetTemporaryNewWeights(std::vector<double> scaling, int first, int last, int off, int i, int j) const {
    std::vector<double> temp(last + 2 - off, 0);
    temp[0] = physical_space_->GetWeights()[first - 1];
    temp[last + 1 - off] = physical_space_->GetWeights()[last + 1];
    while (j - i > 0) {
      double alfi = scaling[i - first];
      double alfj = scaling[j - first];
      temp[i - off] = (physical_space_->GetWeights()[i] - (1 - alfi) * temp[i - off - 1]) / alfi;
      temp[j - off] = (physical_space_->GetWeights()[j] - alfj * temp[j - off + 1]) / (1 - alfj);
      ++i, --j;
    }
    return temp;
  }

  bool IsKnotRemovable(double alfi, std::vector<double> temp, std::vector<double> temp_w, double tolerance,
                       int i, int j, int off, int dimension) const {
    std::array<int, DIM> point_handler_length = physical_space_->GetNumberOfPointsInEachDirection();
    point_handler_length[dimension] = 0;
    util::MultiIndexHandler<DIM> point_handler(point_handler_length);
    int new_control_points =
        physical_space_->GetNumberOfControlPoints() / physical_space_->GetNumberOfPointsInEachDirection()[dimension];
    size_t temp_length = temp.size() / new_control_points;
    for (int l = 0; l < new_control_points; ++l, ++point_handler) {
      size_t offset = l * temp_length;
      std::vector<double> temp1(GetDimension() + 1, temp_w[offset + i - off - 1]);
      std::vector<double> temp2(GetDimension() + 1, temp_w[offset + j - off + 1]);
      for (int k = 0; k < GetDimension(); ++k) {
        temp1[k] = temp[offset + (i - off - 1) * GetDimension() + k] * temp_w[offset + i - off - 1];
        temp2[k] = temp[offset + (j - off + 1) * GetDimension() + k] * temp_w[offset + j - off + 1];
      }
      if (util::VectorUtils<double>::ComputeDistance(temp1, temp2) > tolerance) {
        for (int k = 0; k < GetDimension(); ++k) {
          auto indices = point_handler.GetIndices();
          indices[dimension] = i;
          temp1[k] = physical_space_->GetHomogenousControlPoint(indices).GetValue(k);
          temp2[k] = alfi * temp[offset + (i - off + 1) * GetDimension() + k] * temp_w[offset + i - off + 1]
              + (1 - alfi) * temp[offset + (i - off - 1) * GetDimension() + k] * temp_w[offset + i - off + 1];
        }
        if (util::VectorUtils<double>::ComputeDistance(temp1, temp2) > tolerance) {
//          return false;
        }
      }
    }
    return true;
  }

  std::shared_ptr<WeightedPhysicalSpace<DIM>> physical_space_;
};
}  // namespace spl

#endif  // SRC_SPL_NURBS_H_
