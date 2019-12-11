/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#ifndef SRC_SPL_NURBS_H_
#define SRC_SPL_NURBS_H_

#include <algorithm>
#include <array>
#include <functional>
#include <utility>
#include <vector>

#include "src/spl/b_spline.h"
#include "src/spl/spline.h"
#include "src/spl/nurbs_generator.h"
#include "src/spl/weighted_physical_space.h"

namespace splinelib::src::spl {
template<int PARAMETRIC_DIMENSIONALITY>
class NURBS : public Spline<PARAMETRIC_DIMENSIONALITY> {
 public:
  NURBS(baf::KnotVectors<PARAMETRIC_DIMENSIONALITY> knot_vector, std::array<Degree, PARAMETRIC_DIMENSIONALITY> degree,
        const std::vector<spl::ControlPoint> &control_points,
        std::vector<double> weights) : Spline<PARAMETRIC_DIMENSIONALITY>(knot_vector, degree) {
    std::array<int, PARAMETRIC_DIMENSIONALITY> number_of_points{};
    for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; ++i) {
      number_of_points[i] = knot_vector[i]->GetNumberOfKnots() - degree[i].Get() - 1;
    }
    physical_space_ = std::make_shared<WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>>(
        WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>(control_points, weights, number_of_points));
  }

  explicit NURBS(NURBSGenerator<PARAMETRIC_DIMENSIONALITY> nurbs_generator) : Spline<PARAMETRIC_DIMENSIONALITY>(
      nurbs_generator.GetParameterSpace()) {
    physical_space_ = nurbs_generator.GetWeightedPhysicalSpace();
  }

  NURBS(const NURBS<PARAMETRIC_DIMENSIONALITY> &nurbs) : Spline<PARAMETRIC_DIMENSIONALITY>(nurbs) {
    WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY> weighted_physical_space(*nurbs.physical_space_);
    physical_space_ = std::make_shared<WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>>(weighted_physical_space);
  }

  NURBS(const NURBS<PARAMETRIC_DIMENSIONALITY> &nurbs, const std::vector<spl::ControlPoint> &control_points)
      : Spline<PARAMETRIC_DIMENSIONALITY>(nurbs) {
    WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>
        weighted_physical_space(control_points, nurbs.physical_space_->GetWeights(), nurbs.GetPointsPerDirection());
    physical_space_ = std::make_shared<WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>>(weighted_physical_space);
  }

  NURBS(NURBS<PARAMETRIC_DIMENSIONALITY> &&other) = delete;
  NURBS & operator=(const NURBS<PARAMETRIC_DIMENSIONALITY> &rhs) = delete;
  NURBS & operator=(NURBS<PARAMETRIC_DIMENSIONALITY> &&rhs) = delete;
  ~NURBS() override = default;

  bool AreEqual(const NURBS<PARAMETRIC_DIMENSIONALITY> &rhs, double tolerance =
      util::numeric_settings::GetEpsilon<double>()) const {
    return this->parameter_space_->AreEqual(*rhs.parameter_space_.get(), tolerance)
        && physical_space_->AreEqual(*rhs.physical_space_.get(), tolerance);
  }

  double GetHomogeneousControlPoint(std::array<int, PARAMETRIC_DIMENSIONALITY> indices, int dimension) const {
    return physical_space_->GetHomogenousControlPoint(indices)[Dimension{dimension}];
  }

  void AdjustControlPoints(std::vector<double> scaling, int first, int last, int dimension) override {
    std::array<int, PARAMETRIC_DIMENSIONALITY> point_handler_length = this->GetPointsPerDirection();
    ++point_handler_length[dimension];
    util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> point_handler(point_handler_length);
    std::array<int, PARAMETRIC_DIMENSIONALITY> maximum_point_index =
        physical_space_->GetMaximumPointIndexInEachDirection();
    ++maximum_point_index[dimension];
    point_handler.SetCurrentIndex(maximum_point_index);
    int new_points = physical_space_->GetNumberOfControlPoints() / maximum_point_index[dimension];
    physical_space_->AddControlPoints(new_points);
    for (int i = point_handler.GetNumberOfTotalMultiIndices() - 1; i >= 0; --i, --point_handler) {
      auto current_point = point_handler[Dimension{dimension}];
      std::array<int, PARAMETRIC_DIMENSIONALITY> indices = point_handler.GetCurrentIndex();
      spl::ControlPoint new_control_point = GetNewControlPoint(indices, dimension, scaling, current_point, first, last);
      double new_weight = GetNewWeight(indices, dimension, scaling, current_point, first, last);
      physical_space_->SetControlPoint(indices, new_control_point, dimension,
                                       util::numeric_operations::increment<int>);
      physical_space_->SetWeight(indices, new_weight, dimension, util::numeric_operations::increment<int>);
    }
    physical_space_->IncrementNumberOfPoints(dimension);
  }

  bool RemoveControlPoints(std::vector<double> scaling, int first, int last, int dimension, double tolerance) override {
    int off = first - 1, i = first, j = last;
    std::vector<double> temp_w = GetTempNewWeights(scaling, off, last, i, j, dimension);
    std::vector<double> temp = GetTempNewControlPoints(scaling, temp_w, off, last, i, j, dimension);
    auto diff = static_cast<int>(ceil((j - i) / 2.0));
    i += diff, j -= diff;
    if (!IsKnotRemovable(scaling[i - off - 1], temp, temp_w, tolerance, i, j, off, dimension)) {
      return false;
    }
    SetNewControlPoints(temp, last, i - off, off, dimension);
    SetNewWeights(temp_w, last, i - off, off, dimension);
    physical_space_->RemoveControlPoints(this->GetNumberOfControlPoints() / this->GetPointsPerDirection()[dimension]);
    physical_space_->DecrementNumberOfPoints(dimension);
    return true;
  }

  std::array<std::shared_ptr<spl::NURBS<PARAMETRIC_DIMENSIONALITY>>, 2> SudivideSpline(ParametricCoordinate param_coord,
      int dimension) {
    this->InsertKnot(param_coord, dimension,
                     this->GetDegree(dimension).Get() + 1
                         - this->GetKnotVector(dimension)->GetMultiplicity(param_coord).Get());
    std::array<baf::KnotVectors<PARAMETRIC_DIMENSIONALITY>, 2>
        new_knot_vectors = this->parameter_space_->GetDividedKnotVectors(param_coord, dimension);
    std::array<Degree, PARAMETRIC_DIMENSIONALITY> degrees{};
    for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; ++i) {
      degrees[i] = this->GetDegree(i);
    }
    std::array<std::shared_ptr<spl::NURBS<PARAMETRIC_DIMENSIONALITY>>, 2> subdivided_splines;
    int first = 0;
    for (int i = 0; i < 2; ++i) {
      int length = new_knot_vectors[i][dimension]->GetNumberOfKnots() - degrees[dimension].Get() - 1;
      std::vector<spl::ControlPoint> points = physical_space_->GetDividedControlPoints(first, length, dimension);
      std::vector<double> weights = physical_space_->GetDividedWeights(first, length, dimension);
      spl::NURBS<PARAMETRIC_DIMENSIONALITY> spline(new_knot_vectors[i], degrees, points, weights);
      subdivided_splines[i] = std::make_shared<spl::NURBS<PARAMETRIC_DIMENSIONALITY>>(spline);
      first = length;
    }
    return subdivided_splines;
  }

 private:
  std::shared_ptr<spl::PhysicalSpace<PARAMETRIC_DIMENSIONALITY>> GetPhysicalSpace() const override {
    return physical_space_;
  }

  double GetEvaluatedControlPoint(std::array<ParametricCoordinate, PARAMETRIC_DIMENSIONALITY> param_coord,
                                  std::array<int, PARAMETRIC_DIMENSIONALITY> indices,
                                  int dimension) const override {
    if (this->GetPointDim() == dimension) {
      return this->parameter_space_->GetBasisFunctions(indices, param_coord)
          * physical_space_->GetHomogenousControlPoint(indices)[Dimension{dimension}];
    }
    return this->parameter_space_->GetBasisFunctions(indices, param_coord)
        * physical_space_->GetHomogenousControlPoint(indices)[Dimension{dimension}]
        / GetEvaluatedWeightSum(param_coord);
  }

  double GetEvaluatedDerivativeControlPoint(std::array<ParametricCoordinate, PARAMETRIC_DIMENSIONALITY> param_coord,
                                            std::array<int, PARAMETRIC_DIMENSIONALITY> derivative,
                                            std::array<int, PARAMETRIC_DIMENSIONALITY> indices,
                                            int dimension) const override {
    return GetRationalBasisFunctionDerivative(param_coord, derivative, indices, dimension)
        * physical_space_->GetControlPoint(indices)[Dimension{dimension}];
  }

  double GetRationalBasisFunctionDerivative(std::array<ParametricCoordinate, PARAMETRIC_DIMENSIONALITY> param_coord,
                                            std::array<int, PARAMETRIC_DIMENSIONALITY> derivative,
                                            std::array<int, PARAMETRIC_DIMENSIONALITY> indices,
                                            int dimension) const {
    if (derivative == std::array<int, PARAMETRIC_DIMENSIONALITY>{0}) {
      return this->parameter_space_->GetBasisFunctions(indices, param_coord)
          * physical_space_->GetWeight(indices) / GetEvaluatedWeightSum(param_coord);
    }
    return (GetEvaluatedDerivativeWeight(param_coord, derivative, indices)
        - GetDerivativesSum(param_coord, derivative, indices, dimension))
        / GetEvaluatedWeightSum(param_coord);
  }

  double GetDerivativesSum(std::array<ParametricCoordinate, PARAMETRIC_DIMENSIONALITY> param_coord,
                           std::array<int, PARAMETRIC_DIMENSIONALITY> derivative,
                           std::array<int, PARAMETRIC_DIMENSIONALITY> indices,
                           int dimension) const {
    util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> derivativeHandler(GetDerivativeHandler(derivative));
    derivativeHandler++;
    double sum = 0;
    for (int i = 1; i < derivativeHandler.GetNumberOfTotalMultiIndices(); i++, derivativeHandler++) {
      sum += binomialCoefficient(derivative, derivativeHandler.GetCurrentIndex())
          * GetRationalBasisFunctionDerivative(param_coord,
                                               derivativeHandler.GetComplementaryIndex(),
                                               indices,
                                               dimension)
          * GetEvaluatedDerivativeWeightSum(param_coord, derivativeHandler.GetCurrentIndex());
    }
    return sum;
  }

  double GetEvaluatedWeightSum(std::array<ParametricCoordinate, PARAMETRIC_DIMENSIONALITY> param_coord) const {
    auto first_non_zero = this->GetArrayOfFirstNonZeroBasisFunctions(param_coord);
    util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> basisFunctionHandler(
        this->GetNumberOfBasisFunctionsToEvaluate());
    double sum = 0;
    for (int i = 0; i < basisFunctionHandler.GetNumberOfTotalMultiIndices(); ++i, basisFunctionHandler++) {
      auto indices = basisFunctionHandler.GetCurrentIndex();
      std::transform(indices.begin(), indices.end(), first_non_zero.begin(), indices.begin(), std::plus<>());
      sum += GetEvaluatedWeight(param_coord, indices);
    }
    return sum;
  }

  double GetEvaluatedDerivativeWeightSum(std::array<ParametricCoordinate, PARAMETRIC_DIMENSIONALITY> param_coord,
                                         std::array<int, PARAMETRIC_DIMENSIONALITY> derivative) const {
    auto first_non_zero = this->GetArrayOfFirstNonZeroBasisFunctions(param_coord);
    util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> basisFunctionHandler(
        this->GetNumberOfBasisFunctionsToEvaluate());
    double sum = 0;
    for (int i = 0; i < basisFunctionHandler.GetNumberOfTotalMultiIndices(); ++i, basisFunctionHandler++) {
      auto indices = basisFunctionHandler.GetCurrentIndex();
      std::transform(indices.begin(), indices.end(), first_non_zero.begin(), indices.begin(), std::plus<>());
      sum += GetEvaluatedDerivativeWeight(param_coord, derivative, indices);
    }
    return sum;
  }

  double GetEvaluatedWeight(std::array<ParametricCoordinate, PARAMETRIC_DIMENSIONALITY> param_coord,
      std::array<int, PARAMETRIC_DIMENSIONALITY> indices) const {
    return this->parameter_space_->GetBasisFunctions(indices, param_coord) * physical_space_->GetWeight(indices);
  }

  double GetEvaluatedDerivativeWeight(std::array<ParametricCoordinate, PARAMETRIC_DIMENSIONALITY> param_coord,
                                      std::array<int, PARAMETRIC_DIMENSIONALITY> derivative,
                                      std::array<int, PARAMETRIC_DIMENSIONALITY> indices) const {
    return this->parameter_space_->GetBasisFunctionDerivatives(indices, param_coord, derivative)
        * physical_space_->GetWeight(indices);
  }

  util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> GetDerivativeHandler(
      const std::array<int, PARAMETRIC_DIMENSIONALITY> &derivative) const {
    std::array<int, PARAMETRIC_DIMENSIONALITY> derivative_length{};
    for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; ++i) {
      derivative_length[i] = derivative[i] + 1;
    }
    return util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>(derivative_length);
  }

  int binomialCoefficient(int number, int subset) const {
    if (subset == 0 || subset == number)
      return 1;
    return binomialCoefficient(number - 1, subset - 1) + binomialCoefficient(number - 1, subset);
  }

  int binomialCoefficient(std::array<int, PARAMETRIC_DIMENSIONALITY> numbers,
      std::array<int, PARAMETRIC_DIMENSIONALITY> subsets) const {
    int bc = 1;
    for (int i = 0; i < PARAMETRIC_DIMENSIONALITY; ++i) {
      bc *= binomialCoefficient(numbers[i], subsets[i]);
    }
    return bc;
  }

  spl::ControlPoint GetNewControlPoint(std::array<int, PARAMETRIC_DIMENSIONALITY> indices, int dimension,
      std::vector<double> scaling, int current_point, int first, int last) {
    if (current_point > last) {
      --indices[dimension];
      return physical_space_->GetControlPoint(indices);
    }
    if (current_point >= first) {
      std::array<int, PARAMETRIC_DIMENSIONALITY> lower_indices = indices;
      --lower_indices[dimension];
      spl::ControlPoint upper_control_point = physical_space_->GetHomogenousControlPoint(indices);
      spl::ControlPoint lower_control_point = physical_space_->GetHomogenousControlPoint(lower_indices);
      std::vector<double> coordinates;
      coordinates.reserve(upper_control_point.GetDimensionality());
      double new_weight = GetNewWeight(indices, dimension, scaling, current_point, first, last);
      for (int j = 0; j < upper_control_point.GetDimensionality(); ++j) {
        coordinates.emplace_back((scaling[current_point - first] *
        upper_control_point[Dimension{j}] + (1 - scaling[current_point - first]) *
        lower_control_point[Dimension{j}]) / new_weight);
      }
      return spl::ControlPoint(coordinates);
    }
    return physical_space_->GetControlPoint(indices);
  }

  double GetNewWeight(std::array<int, PARAMETRIC_DIMENSIONALITY> indices, int dimension, std::vector<double> scaling,
                      int current_point, int first, int last) {
    if (current_point > last) {
      --indices[dimension];
      return physical_space_->GetWeight(indices);
    }
    if (current_point >= first) {
      std::array<int, PARAMETRIC_DIMENSIONALITY> lower_indices = indices;
      --lower_indices[dimension];
      double upper_weight = physical_space_->GetWeight(indices);
      double lower_weight = physical_space_->GetWeight(lower_indices);
      return scaling[current_point - first] * upper_weight + (1 - scaling[current_point - first]) * lower_weight;
    }
    return physical_space_->GetWeight(indices);
  }

  void SetNewControlPoints(const std::vector<double> &temp, int last, int ii, int off, int dimension) {
    std::array<int, PARAMETRIC_DIMENSIONALITY> point_handler_length = this->GetPointsPerDirection();
    util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> point_handler(point_handler_length);
    for (int m = 0; m < point_handler.GetNumberOfTotalMultiIndices(); ++m, ++point_handler) {
      int k = point_handler[Dimension{dimension}];
      if (k - off >= 1 && k - off != ii && k < last + 2) {
        int index =
            (point_handler.CollapseDimension(Dimension{dimension}) * (last - off + 2) + k - off) * this->GetPointDim();
        std::vector<double> coordinates(temp.begin() + index, temp.begin() + index + this->GetPointDim());
        auto indices = point_handler.GetCurrentIndex();
        indices[dimension] = k - off < ii ? k : k - 1;
        physical_space_->SetControlPoint(indices, spl::ControlPoint(coordinates), dimension,
                                         util::numeric_operations::decrement<int>);
      }
      if ((k <= off && k - off < 1) || (k >= last + 1 && k < this->GetPointsPerDirection()[dimension])) {
        auto indices = point_handler.GetCurrentIndex();
        indices[dimension] = k <= off ? k : k - 1;
        physical_space_->SetControlPoint(indices, this->GetControlPoint(point_handler.GetCurrentIndex()), dimension,
                                         util::numeric_operations::decrement<int>);
      }
    }
  }

  void SetNewWeights(const std::vector<double> &temp, int last, int ii, int off, int dimension) {
    std::array<int, PARAMETRIC_DIMENSIONALITY> point_handler_length = this->GetPointsPerDirection();
    util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> point_handler(point_handler_length);
    for (int m = 0; m < point_handler.GetNumberOfTotalMultiIndices(); ++m, ++point_handler) {
      int k = point_handler[Dimension{dimension}];
      if (k - off >= 1 && k - off != ii && k < last + 2) {
        int index = point_handler.CollapseDimension(Dimension{dimension}) * (last - off + 2) + k - off;
        auto indices = point_handler.GetCurrentIndex();
        indices[dimension] = k - off < ii ? k : k - 1;
        physical_space_->SetWeight(indices, temp[index], dimension, util::numeric_operations::decrement<int>);
      }
      if ((k <= off && k - off < 1) || (k >= last + 1 && k < this->GetPointsPerDirection()[dimension])) {
        auto indices = point_handler.GetCurrentIndex();
        indices[dimension] = k <= off ? k : k - 1;
        physical_space_->SetWeight(indices, this->GetWeight(point_handler.GetCurrentIndex()), dimension,
                                   util::numeric_operations::decrement<int>);
      }
    }
  }

  std::vector<double> GetTempNewControlPoints(const std::vector<double> &scaling, const std::vector<double> &temp_w,
                                              int off, int last, int i, int j, int dimension) const {
    std::array<int, PARAMETRIC_DIMENSIONALITY> point_handler_length = this->GetPointsPerDirection();
    util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> point_handler(point_handler_length);
    int new_points = this->GetNumberOfControlPoints() / this->GetPointsPerDirection()[dimension];
    std::vector<double> temp(new_points * this->GetPointDim() * (last - off + 2), 0);
    std::shared_ptr<std::vector<double>> temp_ptr = std::make_shared<std::vector<double>>(temp);
    for (int l = 0; l < point_handler.GetNumberOfTotalMultiIndices(); ++l, ++point_handler) {
      if (point_handler[Dimension{dimension}] == off || point_handler[Dimension{dimension}] == last + 1) {
        int diff = point_handler[Dimension{dimension}] == off ? off : last + 1;
        SetTempNewControlPoint(point_handler, temp_ptr, temp_w, 1, diff, off, last, dimension, 0);
      }
    }
    for (; j - i > 0; ++i, --j) {
      point_handler.SetCurrentIndex({0});
      for (int l = 0; l < point_handler.GetNumberOfTotalMultiIndices(); ++l, ++point_handler) {
        if (point_handler[Dimension{dimension}] == i) {
          SetTempNewControlPoint(point_handler, temp_ptr, temp_w, scaling[i - off - 1], i, off, last, dimension, -1);
        }
        if (point_handler[Dimension{dimension}] == j) {
          SetTempNewControlPoint(point_handler, temp_ptr, temp_w, 1 - scaling[j - off - 1], j, off, last, dimension, 1);
        }
      }
    }
    return *temp_ptr;
  }

  void SetTempNewControlPoint(const util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> &point_handler,
                              const std::shared_ptr<std::vector<double>> &temp_ptr, const std::vector<double> &temp_w,
                              double alpha, int x, int off, int last, int dimension, int shift) const {
    int index = point_handler.CollapseDimension(Dimension{dimension}) * (last - off + 2);
    for (int k = 0; k < this->GetPointDim(); ++k) {
      (*temp_ptr)[(index + x - off) * this->GetPointDim() + k] =
          (GetHomogeneousControlPoint(point_handler.GetCurrentIndex(), k)
              - (1 - alpha) * (*temp_ptr)[(index + x - off + shift) * this->GetPointDim() + k]
                  * temp_w[index + x - off + shift]) / alpha / temp_w[index + x - off];
    }
  }

  std::vector<double> GetTempNewWeights(const std::vector<double> &scaling, int off, int last,
                                        int i, int j, int dimension) const {
    std::array<int, PARAMETRIC_DIMENSIONALITY> point_handler_length = this->GetPointsPerDirection();
    util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> point_handler(point_handler_length);
    int new_control_points = this->GetNumberOfControlPoints() / this->GetPointsPerDirection()[dimension];
    std::vector<double> temp_w(new_control_points * (last - off + 2), 0);
    std::shared_ptr<std::vector<double>> temp_w_ptr = std::make_shared<std::vector<double>>(temp_w);
    for (int l = 0; l < point_handler.GetNumberOfTotalMultiIndices(); ++l, ++point_handler) {
      if (point_handler[Dimension{dimension}] == off || point_handler[Dimension{dimension}] == last + 1) {
        int diff = point_handler[Dimension{dimension}] == off ? off : last + 1;
        SetTempNewWeight(point_handler, temp_w_ptr, 1, diff, off, last, dimension, 0);
      }
    }
    for (; j - i > 0; ++i, --j) {
      point_handler.SetCurrentIndex({0});
      for (int l = 0; l < point_handler.GetNumberOfTotalMultiIndices(); ++l, ++point_handler) {
        if (point_handler[Dimension{dimension}] == i) {
          SetTempNewWeight(point_handler, temp_w_ptr, scaling[i - off - 1], i, off, last, dimension, -1);
        }
        if (point_handler[Dimension{dimension}] == j) {
          SetTempNewWeight(point_handler, temp_w_ptr, 1 - scaling[j - off - 1], j, off, last, dimension, 1);
        }
      }
    }
    return *temp_w_ptr;
  }

  void SetTempNewWeight(const util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> &point_handler,
                        const std::shared_ptr<std::vector<double>> &temp_w_ptr,
                        double alpha, int x, int off, int last, int dimension, int shift) const {
    int index = point_handler.CollapseDimension(Dimension{dimension}) * (last - off + 2);
    (*temp_w_ptr)[index + x - off] = (this->GetWeight(point_handler.GetCurrentIndex())
        - (1 - alpha) * (*temp_w_ptr)[index + x - off + shift]) / alpha;
  }

  bool IsKnotRemovable(double alfi, const std::vector<double> &temp, const std::vector<double> &temp_w,
                       double tolerance, int i, int j, int off, int dimension) const {
    auto maxdist = physical_space_->GetMaximumDistanceFromOrigin();
    auto minw = physical_space_->GetMinimumWeight();
    tolerance = tolerance * minw / (1 + maxdist);

    std::array<int, PARAMETRIC_DIMENSIONALITY> point_handler_length = this->GetPointsPerDirection();
    point_handler_length[dimension] = 0;
    util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> point_handler(point_handler_length);
    int new_points = this->GetNumberOfControlPoints() / this->GetPointsPerDirection()[dimension];
    size_t temp_length = temp.size() / new_points;
    size_t temp_w_length = temp_w.size() / new_points;
    for (int l = 0; l < new_points; ++l, ++point_handler) {
      size_t offset = l * temp_length;
      size_t w_offset = l * temp_w_length;
      std::vector<double> temp1(this->GetPointDim() + 1, temp_w[w_offset + i - off - 1]);
      std::vector<double> temp2(this->GetPointDim() + 1, temp_w[w_offset + j - off + 1]);
      for (int k = 0; k < this->GetPointDim(); ++k) {
        temp1[k] = temp[offset + (i - off - 1) * this->GetPointDim() + k] * temp_w[w_offset + i - off - 1];
        temp2[k] = temp[offset + (j - off + 1) * this->GetPointDim() + k] * temp_w[w_offset + j - off + 1];
      }
      if (util::vector_utils::ComputeDistance(temp1, temp2) > tolerance) {
        auto indices = point_handler.GetCurrentIndex();
        indices[dimension] = i;
        for (int k = 0; k < this->GetPointDim(); ++k) {
          temp1[k] = physical_space_->GetHomogenousControlPoint(indices)[Dimension{k}];
          temp2[k] = alfi * temp[offset + (i - off + 1) * this->GetPointDim() + k] * temp_w[w_offset + i - off + 1]
              + (1 - alfi) * temp[offset + (i - off - 1) * this->GetPointDim() + k] * temp_w[w_offset + i - off - 1];
        }
        temp1[this->GetPointDim()] = physical_space_->GetWeight(indices);
        temp2[this->GetPointDim()] =
            alfi * temp_w[w_offset + i - off + 1] + (1 - alfi) * temp_w[w_offset + i - off - 1];
        if (util::vector_utils::ComputeDistance(temp1, temp2) > tolerance) {
          return false;
        }
      }
    }
    return true;
  }

  void SetNewControlPoint(spl::ControlPoint control_point, double weight,
      std::array<int, PARAMETRIC_DIMENSIONALITY> indices) override {
    physical_space_->SetWeightedControlPoint(indices, control_point, weight);
  }

  std::shared_ptr<WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>> physical_space_;
};
}  // namespace splinelib::src::spl

#endif  // SRC_SPL_NURBS_H_
