/* Copyright 2019 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

This file is part of SplineLib.

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
<http://www.gnu.org/licenses/>.*/

#include <vector>


template<int PARAMETRIC_DIMENSIONALITY>
WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>::WeightedPhysicalSpace(
    std::vector<spl::ControlPoint> const &control_points, std::vector<double> const &weights,
    std::array<int, PARAMETRIC_DIMENSIONALITY> const &number_of_points)
    : PhysicalSpace<PARAMETRIC_DIMENSIONALITY>(control_points, number_of_points), weights_(std::move(weights)) {
  if (control_points.size() != weights_.size()) {
    throw std::runtime_error("The number of control points and weights has to be the same.");
  }
}

template<int PARAMETRIC_DIMENSIONALITY>
WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>::WeightedPhysicalSpace(WeightedPhysicalSpace const &other)
    : PhysicalSpace<PARAMETRIC_DIMENSIONALITY>(other) {
  this->weights_ = other.weights_;
}

template<int PARAMETRIC_DIMENSIONALITY>
bool WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>::AreEqual(WeightedPhysicalSpace const &rhs, double tolerance) const {
  return std::equal(weights_.begin(), weights_.end(),
                    rhs.weights_.begin(), rhs.weights_.end(),
                    [&](double weight_a, double weight_b) {
                      return util::numeric_settings::AreEqual<double>(weight_a, weight_b, tolerance);
                    }) && PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::AreEqual(rhs, Tolerance{tolerance});
}

template<int PARAMETRIC_DIMENSIONALITY>
ControlPoint WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>::GetHomogenousControlPoint(
    std::array<int, PARAMETRIC_DIMENSIONALITY> const &multi_index) const {
  auto point_handler = util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>(this->number_of_points_per_dimension_);
  point_handler.SetCurrentIndex(multi_index);
  int first = this->dimensionality_ * point_handler.GetCurrent1DIndex();
  std::vector<double> coordinates;
  coordinates.reserve(this->dimensionality_);
  for (int coordinate = 0; coordinate < this->dimensionality_; coordinate++) {
    coordinates.emplace_back(this->control_points_[first + coordinate] * weights_[first / this->dimensionality_]);
  }
  coordinates.push_back(this->weights_[point_handler.GetCurrent1DIndex()]);
  return spl::ControlPoint(coordinates);
}

template<int PARAMETRIC_DIMENSIONALITY>
ControlPoint WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>::GetHomogenousControlPoint(int index_1d) const {
  auto point_handler = util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>(this->number_of_points_per_dimension_);
  point_handler.SetCurrentIndex(index_1d);
  int first = this->dimensionality_ * point_handler.GetCurrent1DIndex();
  std::vector<double> coordinates;
  coordinates.reserve(this->dimensionality_);
  for (int coordinate = 0; coordinate < this->dimensionality_; coordinate++) {
    coordinates.emplace_back(this->control_points_[first + coordinate] * weights_[first / this->dimensionality_]);
  }
  coordinates.push_back(this->weights_[point_handler.GetCurrent1DIndex()]);
  return spl::ControlPoint(coordinates);
}

template<int PARAMETRIC_DIMENSIONALITY>
Weight WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>::GetWeight(
    std::array<int, PARAMETRIC_DIMENSIONALITY> const &multi_index) const {
  auto point_handler = util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>(this->number_of_points_per_dimension_);
  point_handler.SetCurrentIndex(multi_index);
  int first = point_handler.GetCurrent1DIndex();
  return Weight{weights_[first]};
}

template<int PARAMETRIC_DIMENSIONALITY>
Weight WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>::GetWeight(int index_1d) const {
  auto point_handler = util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>(this->number_of_points_per_dimension_);
  point_handler.SetCurrentIndex(index_1d);
  int first = point_handler.GetCurrent1DIndex();
  return Weight{weights_[first]};
}

template<int PARAMETRIC_DIMENSIONALITY>
double WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>::GetMinimumWeight() const {
  double minimum = weights_[0];
  for (const auto &weight : weights_) {
    if (weight < minimum) minimum = weight;
  }
  return minimum;
}

template<int PARAMETRIC_DIMENSIONALITY>
void WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>::SetWeightedControlPoint(
    std::array<int, PARAMETRIC_DIMENSIONALITY> const &multi_index,
    ControlPoint const &control_point, double weight) {
  PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::SetControlPoint(multi_index, control_point);
  SetWeight(multi_index, weight);
}

template<int PARAMETRIC_DIMENSIONALITY>
void WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>::SetWeight(
    std::array<int, PARAMETRIC_DIMENSIONALITY> const &multi_index, double weight, int dimension, int (*before)(int)) {
  const std::array<int, PARAMETRIC_DIMENSIONALITY> number_of_points_before(this->number_of_points_per_dimension_);
  if (before != nullptr)
    this->number_of_points_per_dimension_[dimension] = before(this->number_of_points_per_dimension_[dimension]);
  auto point_handler = util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>(this->number_of_points_per_dimension_);
  point_handler.SetCurrentIndex(multi_index);
  int first = point_handler.GetCurrent1DIndex();
  weights_[first] = weight;
  this->number_of_points_per_dimension_ = number_of_points_before;
}

template<int PARAMETRIC_DIMENSIONALITY>
void WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>::SetWeight(int index_1d, double weight, int dimension,
               int (*before)(int)) {
  const std::array<int, PARAMETRIC_DIMENSIONALITY> number_of_points_before(this->number_of_points_per_dimension_);
  if (before != nullptr)
    this->number_of_points_per_dimension_[dimension] = before(this->number_of_points_per_dimension_[dimension]);
  auto point_handler = util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY>(this->number_of_points_per_dimension_);
  point_handler.SetCurrentIndex(index_1d);
  int first = point_handler.GetCurrent1DIndex();
  weights_[first] = weight;
  this->number_of_points_per_dimension_ = number_of_points_before;
}

template<int PARAMETRIC_DIMENSIONALITY>
void WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>::AddControlPoints(int number) {
  PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::AddControlPoints(number);
  for (int i = 0; i < number; ++i) {
    weights_.emplace_back(0.0);
  }
}

template<int PARAMETRIC_DIMENSIONALITY>
void WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>::RemoveControlPoints(int number) {
  PhysicalSpace<PARAMETRIC_DIMENSIONALITY>::RemoveControlPoints(number);
  weights_.erase(weights_.end() - number, weights_.end());
}

template<int PARAMETRIC_DIMENSIONALITY>
std::vector<double> WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>::GetWeights() const {
  return weights_;
}

template<int PARAMETRIC_DIMENSIONALITY>
std::vector<double> WeightedPhysicalSpace<PARAMETRIC_DIMENSIONALITY>::GetDividedWeights(
    int first, int length, int dimension) {
  std::vector<double> weights;
  std::array<int, PARAMETRIC_DIMENSIONALITY> point_handler_length = this->GetNumberOfPointsPerDirection();
  point_handler_length[dimension] = length;
  util::MultiIndexHandler<PARAMETRIC_DIMENSIONALITY> point_handler(point_handler_length);
  for (int i = 0; i < point_handler.GetNumberOfTotalMultiIndices(); ++i, ++point_handler) {
    auto indices = point_handler.GetCurrentIndex();
    indices[dimension] += first;
    weights.emplace_back(GetWeight(indices).Get());
  }
  return weights;
}
